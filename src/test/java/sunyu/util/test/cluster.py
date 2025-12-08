import dash
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc, html, Output, Input
from sklearn.cluster import DBSCAN

# 读取实际数据文件
try:
    # 读取高斯投影坐标数据
    data_path = "D:/GitLab/util-gis/testFiles/gauss_xy.txt"
    coordinates = []
    with open(data_path, 'r', encoding='utf-8') as f:
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
        df = pd.DataFrame(coordinates, columns=['x', 'y'])
        print(f"成功读取数据文件，共 {len(df)} 个点")
    else:
        print("数据文件为空或格式错误，使用模拟数据")
        # 回退到模拟数据
        np.random.seed(42)
        line1 = np.column_stack([np.linspace(0, 100, 50), np.random.normal(0, 0.1, 50) + 0])
        line2 = np.column_stack([np.linspace(0, 100, 50), np.random.normal(0, 0.1, 50) + 3])
        road = np.column_stack([np.arange(0, 100, 10), np.arange(0, 100, 10) + 50])
        data = np.vstack([line1, line2, road])
        df = pd.DataFrame(data, columns=['x', 'y'])

except FileNotFoundError:
    print("数据文件未找到，使用模拟数据")
    # 回退到模拟数据
    np.random.seed(42)
    line1 = np.column_stack([np.linspace(0, 100, 50), np.random.normal(0, 0.1, 50) + 0])
    line2 = np.column_stack([np.linspace(0, 100, 50), np.random.normal(0, 0.1, 50) + 3])
    road = np.column_stack([np.arange(0, 100, 10), np.arange(0, 100, 10) + 50])
    data = np.vstack([line1, line2, road])
    df = pd.DataFrame(data, columns=['x', 'y'])

# 打印数据范围信息
print(f"数据范围: X ∈ [{df['x'].min():.2f}, {df['x'].max():.2f}], Y ∈ [{df['y'].min():.2f}, {df['y'].max():.2f}]")

# 初始化Dash应用
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container([
    # 移除标题和数据文件信息展示

    # 参数输入和按钮在同一行 - 紧凑型布局
    dbc.Row([
        dbc.Col([
            html.Label("邻域半径 eps (米):", style={'margin-bottom': '5px'}),
            dbc.Input(id='eps-input', type='number', value=5, min=0.01, max=50, step=0.01,
                      placeholder="输入邻域半径", className="mb-1", size="sm"),
            html.Div(id='eps-value', className="mt-1 small")
        ], width=4),

        dbc.Col([
            html.Label("最小点数 minPts:", style={'margin-bottom': '5px'}),
            dbc.Input(id='minpts-input', type='number', value=3, min=2, max=50, step=1,
                      placeholder="输入最小点数", className="mb-1", size="sm"),
            html.Div(id='minpts-value', className="mt-1 small")
        ], width=4),

        dbc.Col([
            html.Label("\u00A0", style={'margin-bottom': '5px'}),  # 空标签用于对齐
            dbc.Button("执行聚类", id='cluster-button', color="primary", className="mb-1 w-100", size="sm"),
        ], width=4)
    ], className="mb-2"),

    # 聚类结果图 - 全屏显示
    dcc.Graph(id='cluster-graph', style={'height': 'calc(100vh - 120px)', 'width': '100%'}),

    # 结果统计 - 紧凑型
    dbc.Row([
        dbc.Col(html.Div(id='cluster-stats'), width=12)
    ], className="mt-2")
], fluid=True, style={'height': '100vh', 'padding': '10px'})  # fluid=True 让容器占满宽度，设置全屏高度和紧凑内边距


@app.callback(
    [Output('eps-value', 'children'),
     Output('minpts-value', 'children'),
     Output('cluster-graph', 'figure'),
     Output('cluster-stats', 'children')],
    [Input('cluster-button', 'n_clicks')],
    [dash.dependencies.State('eps-input', 'value'),
     dash.dependencies.State('minpts-input', 'value')]
)
def update_clustering(n_clicks, eps, minPts):
    # 如果没有点击按钮，显示所有点为黑色
    if n_clicks is None:
        # 初始状态显示所有点为黑色
        initial_fig = go.Figure()
        initial_fig.add_trace(go.Scatter(
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

        initial_fig.update_layout(
            title='原始数据分布 - 点击"执行聚类"按钮开始分析',
            xaxis_title='X 坐标 (米)',
            yaxis_title='Y 坐标 (米)',
            showlegend=False,
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
        return '', '', initial_fig, None

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
        return f'eps: {eps:.2f}m' if eps is not None else 'eps: 未设置', f'minPts: {minPts}' if minPts is not None else 'minPts: 未设置', empty_fig, empty_stats

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
        title=f'聚类结果 (eps={eps}, minPts={minPts})',
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
        html.P(f"总点数: {len(df)}"),
        html.Hr(),
        html.P(f"当前参数: eps = {eps:.2f}米, minPts = {minPts}", className="mb-0")
    ], color="info")

    return f'eps: {eps:.2f}m', f'minPts: {minPts}', fig, stats


if __name__ == '__main__':
    app.run(debug=True, port=8050)