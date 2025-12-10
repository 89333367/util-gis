import dash
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc, html, Output, Input, State
from sklearn.cluster import DBSCAN
import os

# 全局变量存储数据
df = None

# 文件读取函数
def read_coordinate_file(file_path):
    """读取坐标文件"""
    coordinates = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:  # 跳过空行
                    parts = line.split(',')
                    if len(parts) >= 2:
                        try:
                            # 处理科学计数法
                            x = float(parts[0].strip())
                            y = float(parts[1].strip())
                            coordinates.append([x, y])
                        except ValueError:
                            continue
        
        if coordinates:
            df_data = pd.DataFrame(coordinates, columns=['x', 'y'])
            print(f"成功读取数据文件，共 {len(df_data)} 个点")
            print(f"数据范围: X ∈ [{df_data['x'].min():.2f}, {df_data['x'].max():.2f}], Y ∈ [{df_data['y'].min():.2f}, {df_data['y'].max():.2f}]")
            return df_data
        else:
            print("数据文件为空或格式错误")
            return None
    except Exception as e:
        print(f"读取文件失败: {e}")
        return None

# 初始化Dash应用
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# 添加自定义CSS样式
app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>DBSCAN聚类分析</title>
        {%favicon%}
        {%css%}
        <style>
            .control-table {
                width: 100%;
                margin-bottom: 10px;
                border-collapse: collapse;
            }
            .control-table td {
                padding: 2px;
                vertical-align: middle;
            }
            .control-table input, .control-table button {
                width: 100%;
                box-sizing: border-box;
            }
            .param-display {
                font-size: 10px;
                line-height: 1.2;
                color: #6c757d;
                padding-left: 5px;
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''

app.layout = dbc.Container([
    # 使用HTML表格实现真正的一行布局
    html.Table([
        html.Tr([
            html.Td([
                dbc.Button("选择数据文件", id='select-file-button', color="primary", size="sm")
            ], style={'width': '15%'}),
            html.Td([
                html.Div([
                    html.Label("eps:", style={'fontSize': '10px', 'marginRight': '2px'}),
                    dbc.Input(id='eps-input', type='number', value=5, min=1, max=50, step=0.01,
                             placeholder="eps", size="sm", style={'width': '60px'})
                ], style={'display': 'flex', 'alignItems': 'center'})
            ], style={'width': '12%'}),
            html.Td([
                html.Div([
                    html.Label("minPts:", style={'fontSize': '10px', 'marginRight': '2px'}),
                    dbc.Input(id='minpts-input', type='number', value=20, min=2, max=50, step=1,
                             placeholder="minPts", size="sm", style={'width': '60px'})
                ], style={'display': 'flex', 'alignItems': 'center'})
            ], style={'width': '12%'}),
            html.Td([
                dbc.Button("执行聚类", id='cluster-button', color="primary", size="sm")
            ], style={'width': '16%'})
        ])
    ], className="control-table"),

    # 聚类结果图 - 全屏显示
    dcc.Graph(id='cluster-graph', 
              figure={
                  'data': [],
                  'layout': {
                      'title': '请选择文件并点击"加载数据"按钮',
                      'xaxis': {'title': 'X 坐标 (米)'},
                      'yaxis': {'title': 'Y 坐标 (米)'}
                  }
              },
              style={'height': 'calc(100vh - 120px)', 'width': '100%'}),

    # 结果统计 - 紧凑型
    dbc.Row([
        dbc.Col(html.Div(id='cluster-stats'), width=12)
    ], className="mt-2")
], fluid=True, style={'height': '100vh', 'padding': '10px'})  # fluid=True 让容器占满宽度，设置全屏高度和紧凑内边距


# 文件选择回调（直接加载数据）
@app.callback(
    [Output('cluster-graph', 'figure'),
     Output('cluster-stats', 'children')],
    [Input('select-file-button', 'n_clicks')],
    prevent_initial_call=True
)
def select_and_load_file(n_clicks):
    if n_clicks is None:
        return dash.no_update, dash.no_update
    
    # 使用tkinter打开文件选择对话框
    import tkinter as tk
    from tkinter import filedialog
    
    # 创建隐藏的根窗口
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    
    # 打开文件选择对话框
    file_path = filedialog.askopenfilename(
        title="选择坐标文件",
        filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")]
    )
    
    # 销毁窗口
    root.destroy()
    
    if not file_path:
        return dash.no_update, dash.no_update
    
    # 读取文件
    global df
    df = read_coordinate_file(file_path)
    
    if df is None:
        error_fig = go.Figure()
        error_fig.update_layout(
            title='文件读取失败，请检查文件路径和格式',
            xaxis_title='X 坐标 (米)',
            yaxis_title='Y 坐标 (米)'
        )
        error_stats = dbc.Alert([
            html.H4("文件读取错误", className="alert-heading"),
            html.P("请确保文件存在且格式正确（每行包含x,y坐标，用逗号分隔）")
        ], color="danger")
        return error_fig, error_stats
    
    # 显示黑色数据点
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df['x'],
        y=df['y'],
        mode='markers',
        marker=dict(
            color='black',
            size=8,
            opacity=0.7
        ),
        name='原始数据',
        hovertemplate='<b>原始数据点</b><br>X: %{x:.2f}<br>Y: %{y:.2f}<extra></extra>'
    ))

    fig.update_layout(
        xaxis_title='X 坐标 (米)',
        yaxis_title='Y 坐标 (米)',
        showlegend=False,
        yaxis=dict(scaleanchor="x", scaleratio=1),
        margin=dict(l=50, r=50, t=50, b=50),
        autosize=True
    )
    
    stats = dbc.Alert([
        html.H4("数据加载成功", className="alert-heading"),
        html.P(f"数据点数量: {len(df)}"),
        html.P(f"X坐标范围: [{df['x'].min():.2f}, {df['x'].max():.2f}]"),
        html.P(f"Y坐标范围: [{df['y'].min():.2f}, {df['y'].max():.2f}]")
    ], color="success")
    
    return fig, stats


 
 # 聚类分析回调
@app.callback(
    [Output('cluster-graph', 'figure', allow_duplicate=True),
     Output('cluster-stats', 'children', allow_duplicate=True)],
    [Input('cluster-button', 'n_clicks')],
    [State('eps-input', 'value'),
     State('minpts-input', 'value')],
    prevent_initial_call=True
)
def update_clustering(n_clicks, eps, minPts):
    # 检查是否有数据
    if df is None:
        empty_fig = go.Figure()
        empty_fig.update_layout(
            title='请先选择文件并加载数据',
            xaxis_title='X 坐标 (米)',
            yaxis_title='Y 坐标 (米)'
        )
        return empty_fig, None
    
    # 输入验证
    if eps is None or minPts is None or eps <= 0 or minPts < 2:
        # 返回空图形
        empty_fig = go.Figure()
        empty_fig.update_layout(
            title='请输入有效的参数 (eps > 0, minPts ≥ 2)',
            xaxis_title='X 坐标 (米)',
            yaxis_title='Y 坐标 (米)',
            # 设置等比例坐标轴
            yaxis=dict(
                scaleanchor="x",
                scaleratio=1,
            ),
            # 自适应边距
            margin=dict(l=50, r=50, t=50, b=50),
            # 让图形自适应容器大小
            autosize=True
        )
        empty_stats = dbc.Alert([
            html.H4("参数错误", className="alert-heading"),
            html.P("请确保 eps > 0 且 minPts ≥ 2")
        ], color="warning")
        return empty_fig, empty_stats

    # 输入验证
    if eps is None or minPts is None or eps <= 0 or minPts < 2:
        # 返回空图形
        empty_fig = go.Figure()
        empty_fig.update_layout(
            title='请输入有效的参数 (eps > 0, minPts ≥ 2)',
            xaxis_title='X 坐标 (米)',
            yaxis_title='Y 坐标 (米)',
            # 设置等比例坐标轴
            yaxis=dict(
                scaleanchor="x",
                scaleratio=1,
            ),
            # 自适应边距
            margin=dict(l=50, r=50, t=50, b=50),
            # 让图形自适应容器大小
            autosize=True
        )
        empty_stats = dbc.Alert([
            html.H4("参数错误", className="alert-heading"),
            html.P("请确保 eps > 0 且 minPts ≥ 2")
        ], color="warning")
        return empty_fig, empty_stats

    # 执行DBSCAN
    dbscan = DBSCAN(eps=eps, min_samples=minPts)
    labels = dbscan.fit_predict(df[['x', 'y']])

    # 统计
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)

    # 颜色映射
    unique_labels = set(labels)
    colors = px.colors.qualitative.Set1 + px.colors.qualitative.Set2
    color_map = {label: colors[i % len(colors)] for i, label in enumerate(unique_labels)}
    color_map[-1] = '#000000'  # 噪声黑色

    # 创建图形
    fig = go.Figure()

    for label in unique_labels:
        mask = labels == label
        color = color_map[label]
        cluster_data = df[mask]
        point_count = len(cluster_data)
        name = f'簇 {label} ({point_count}点)' if label != -1 else f'噪声 ({point_count}点)'

        # 绘制聚类点 - 增强显示效果
        marker_size = 10 if label != -1 else 8  # 增大点的大小
        marker_opacity = 0.8 if label != -1 else 0.6  # 噪声点稍微透明

        fig.add_trace(go.Scatter(
            x=cluster_data['x'],
            y=cluster_data['y'],
            mode='markers',
            marker=dict(
                color=color,
                size=marker_size,
                opacity=marker_opacity,
                line=dict(width=1, color='white')  # 添加白色边框增强对比度
            ),
            name=name,
            hovertemplate=f'<b>{name}</b><br>X: %{{x:.2f}}<br>Y: %{{y:.2f}}<br>点数: {point_count}<extra></extra>'
        ))

        # 移除轮廓边界，只使用颜色区分聚类
        # 镂空数据不适合凸壳边界，仅用颜色区分不同类别

    fig.update_layout(
        xaxis_title='X 坐标 (米)',
        yaxis_title='Y 坐标 (米)',
        showlegend=True,
        hovermode='closest',
        # 设置等比例坐标轴
        yaxis=dict(
            scaleanchor="x",
            scaleratio=1,
        ),
        # 自适应边距
        margin=dict(l=50, r=50, t=50, b=50),
        # 让图形自适应容器大小
        autosize=True
    )

    # 统计文本 - 增强信息显示
    stats = dbc.Alert([
        html.H4(f"聚类统计", className="alert-heading"),
        html.P(f"簇数量: {n_clusters}"),
        html.P(f"噪声点: {n_noise} ({n_noise / len(df) * 100:.1f}%)"),
        html.P(f"总点数: {len(df)}")
    ], color="info")

    return fig, stats


if __name__ == '__main__':
    app.run(debug=True, port=8050)