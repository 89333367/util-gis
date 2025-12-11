package sunyu.util;

import cn.hutool.core.io.FileUtil;
import cn.hutool.core.util.StrUtil;
import cn.hutool.log.Log;
import cn.hutool.log.LogFactory;
import elki.clustering.dbscan.DBSCAN;
import elki.data.Cluster;
import elki.data.Clustering;
import elki.data.DoubleVector;
import elki.data.model.Model;
import elki.data.type.TypeUtil;
import elki.database.Database;
import elki.database.StaticArrayDatabase;
import elki.database.ids.DBIDIter;
import elki.database.ids.DBIDRange;
import elki.database.relation.Relation;
import elki.datasource.ArrayAdapterDatabaseConnection;
import elki.distance.minkowski.EuclideanDistance;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.AbstractDocument;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.DocumentFilter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.*;
import java.util.List;

/**
 * DBSCAN聚类可视化工具
 * 支持高斯投影坐标文件加载、参数配置和聚类结果可视化
 *
 * @author SunYu
 */
public class ClusterVisualizationGui extends JFrame {
    private static final Log log = LogFactory.get();

    // GUI组件
    private JTextField filePathField;
    private JTextField epsField;
    private JTextField minPtsField;
    private JButton selectFileButton;
    private JButton clusterButton;
    private ChartPanel chartPanel;
    private JFreeChart chart;
    private JPanel legendPanel;  // 自定义图例面板

    // 数据存储
    private List<double[]> coordinates = new ArrayList<>();
    private XYSeriesCollection dataset;

    // 聚类颜色配置
    private static final Color[] CLUSTER_COLORS = {new Color(255, 99, 71), new Color(30, 144, 255), new Color(50, 205, 50), new Color(255, 215, 0), new Color(138, 43, 226), new Color(255, 140, 0), new Color(220, 20, 60), new Color(0, 191, 255), new Color(255, 20, 147), new Color(0, 250, 154), new Color(255, 69, 0), new Color(106, 90, 205), new Color(255, 255, 0), new Color(186, 85, 211), new Color(0, 255, 127), new Color(123, 104, 238)};
    private static final Color NOISE_COLOR = Color.BLACK;

    /**
     * 构造函数：初始化聚类可视化工具的GUI组件
     */
    public ClusterVisualizationGui() {
        initializeComponents();
        setupLayout();
        setupEventListeners();
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setTitle("空间密集聚类测试工具 - SunYu");
        setSize(1200, 800);
        setLocationRelativeTo(null); // 居中显示
    }

    /**
     * 初始化GUI组件
     * 创建并配置所有界面元素，包括文件选择、参数输入、图表和图例等
     */
    private void initializeComponents() {
        // 文件选择组件
        filePathField = new JTextField(100);
        filePathField.setEditable(false);

        selectFileButton = new JButton("选择数据文件");

        // eps参数输入验证 - 必须大于1.0的浮点数
        epsField = new JTextField("5", 8);
        ((AbstractDocument) epsField.getDocument()).setDocumentFilter(new DocumentFilter() {
            /**
             * 插入字符串时验证eps参数的有效性
             * 只有当新字符串通过验证时才允许插入
             */
            @Override
            public void insertString(FilterBypass fb, int offset, String string, AttributeSet attr) throws BadLocationException {
                if (string == null) return;
                String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
                String newStr = currentText.substring(0, offset) + string + currentText.substring(offset);
                if (isValidEps(newStr)) {
                    super.insertString(fb, offset, string, attr);
                }
            }

            /**
             * 替换字符串时验证eps参数的有效性
             * 只有当新字符串通过验证时才允许替换
             */
            @Override
            public void replace(FilterBypass fb, int offset, int length, String text, AttributeSet attrs) throws BadLocationException {
                if (text == null) return;
                String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
                String newStr = currentText.substring(0, offset) + text + currentText.substring(offset + length);
                if (isValidEps(newStr)) {
                    super.replace(fb, offset, length, text, attrs);
                }
            }

            /**
             * 删除字符串时验证eps参数的有效性
             * 只有当删除后的新字符串通过验证时才允许删除
             */
            @Override
            public void remove(FilterBypass fb, int offset, int length) throws BadLocationException {
                String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
                String newStr = currentText.substring(0, offset) + currentText.substring(offset + length);
                if (isValidEps(newStr)) {
                    super.remove(fb, offset, length);
                }
            }

            /**
             * 验证eps参数的有效性
             * 检查输入文本是否为有效的浮点数且大于等于1.0
             *
             * @param text 要验证的文本
             * @return 如果文本有效返回true，否则返回false
             */
            private boolean isValidEps(String text) {
                if (text.isEmpty()) return true; // 允许空字符串

                // 检查是否只包含数字、小数点和负号
                for (int i = 0; i < text.length(); i++) {
                    char c = text.charAt(i);
                    if (!Character.isDigit(c) && c != '.' && c != '-') {
                        return false;
                    }
                }

                // 检查小数点数量
                int dotCount = 0;
                for (int i = 0; i < text.length(); i++) {
                    if (text.charAt(i) == '.') {
                        dotCount++;
                    }
                }
                if (dotCount > 1) return false;

                // 检查负号位置
                if (text.indexOf('-') > 0) return false;

                try {
                    float value = Float.parseFloat(text);
                    return value >= 1.0f; // eps参数必须大于等于1.0
                } catch (NumberFormatException e) {
                    return false;
                }
            }
        });

        minPtsField = new JTextField("20", 8);
        // minPts参数输入验证 - 必须是正整数
        ((AbstractDocument) minPtsField.getDocument()).setDocumentFilter(new DocumentFilter() {
            /**
             * 插入字符串时验证minPts参数的有效性
             * 只有当新字符串通过验证时才允许插入
             */
            @Override
            public void insertString(FilterBypass fb, int offset, String string, AttributeSet attr) throws BadLocationException {
                if (string == null) return;
                String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
                String newStr = currentText.substring(0, offset) + string + currentText.substring(offset);
                if (isValidMinPts(newStr)) {
                    super.insertString(fb, offset, string, attr);
                }
            }

            /**
             * 替换字符串时验证minPts参数的有效性
             * 只有当新字符串通过验证时才允许替换
             */
            @Override
            public void replace(FilterBypass fb, int offset, int length, String text, AttributeSet attrs) throws BadLocationException {
                if (text == null) return;
                String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
                String newStr = currentText.substring(0, offset) + text + currentText.substring(offset + length);
                if (isValidMinPts(newStr)) {
                    super.replace(fb, offset, length, text, attrs);
                }
            }

            /**
             * 删除字符串时验证minPts参数的有效性
             * 只有当删除后的新字符串通过验证时才允许删除
             */
            @Override
            public void remove(FilterBypass fb, int offset, int length) throws BadLocationException {
                String currentText = fb.getDocument().getText(0, fb.getDocument().getLength());
                String newStr = currentText.substring(0, offset) + currentText.substring(offset + length);
                if (isValidMinPts(newStr)) {
                    super.remove(fb, offset, length);
                }
            }

            /**
             * 验证minPts参数的有效性
             * 检查输入文本是否为有效的正整数且大于等于2
             *
             * @param text 要验证的文本
             * @return 如果文本有效返回true，否则返回false
             */
            private boolean isValidMinPts(String text) {
                if (text.isEmpty()) return true; // 允许空字符串

                // 检查是否只包含数字
                for (int i = 0; i < text.length(); i++) {
                    char c = text.charAt(i);
                    if (!Character.isDigit(c)) {
                        return false;
                    }
                }

                try {
                    int value = Integer.parseInt(text);
                    return value >= 2; // minPts参数必须大于等于2
                } catch (NumberFormatException e) {
                    return false;
                }
            }
        });

        clusterButton = new JButton("执行聚类");
        clusterButton.setEnabled(false);

        // 初始化图表
        dataset = new XYSeriesCollection();
        chart = ChartFactory.createScatterPlot("", "", "", dataset, PlotOrientation.VERTICAL, true, true, false);
        customizeChart();
        chart.removeLegend(); // 使用自定义图例

        // 初始化图例面板
        legendPanel = new JPanel();
        legendPanel.setLayout(new BoxLayout(legendPanel, BoxLayout.Y_AXIS));
        legendPanel.setPreferredSize(new Dimension(200, 600));
        legendPanel.setBorder(BorderFactory.createTitledBorder("图例"));
        legendPanel.setBackground(Color.WHITE);

        legendPanel.removeAll();
        JLabel emptyLabel = new JLabel("暂无数据");
        emptyLabel.setFont(new Font(Font.SANS_SERIF, Font.ITALIC, 12));
        emptyLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        legendPanel.add(Box.createVerticalGlue());
        legendPanel.add(emptyLabel);
        legendPanel.add(Box.createVerticalGlue());

        // 初始化图表面板
        chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        chartPanel.setMouseWheelEnabled(false);
        chartPanel.setMouseZoomable(false);
        chartPanel.setFillZoomRectangle(false);  // 禁用内置缩放矩形，使用自定义框选
        chartPanel.setMaximumDrawWidth(32768);
        chartPanel.setMaximumDrawHeight(32768);
        chartPanel.setMinimumDrawWidth(100);
        chartPanel.setMinimumDrawHeight(100);
        chartPanel.setDomainZoomable(false);    // 禁用内置X轴缩放，使用自定义框选
        chartPanel.setRangeZoomable(false);    // 禁用内置Y轴缩放，使用自定义框选
        chartPanel.setPopupMenu(null);
        chartPanel.setRefreshBuffer(true);
        chartPanel.setDoubleBuffered(true);

        setupUnifiedMouseHandler(chartPanel);  // 使用统一的鼠标处理器
        setupMouseHoverZoom(chartPanel);
    }

    /**
     * 配置图表样式
     * 设置图表的背景、坐标轴、渲染器等视觉效果
     */
    private void customizeChart() {
        XYPlot plot = chart.getXYPlot();

        // 设置背景色
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);

        // 隐藏网格线，让图表更简洁
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);

        // 设置坐标轴
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();

        domainAxis.setAutoRangeIncludesZero(false);
        rangeAxis.setAutoRangeIncludesZero(false);

        // 隐藏坐标轴标签和刻度
        domainAxis.setTickLabelsVisible(false);  // 隐藏X轴刻度标签
        domainAxis.setTickMarksVisible(false);   // 隐藏X轴刻度线
        domainAxis.setAxisLineVisible(false);    // 隐藏X轴线

        rangeAxis.setTickLabelsVisible(false);   // 隐藏Y轴刻度标签
        rangeAxis.setTickMarksVisible(false);    // 隐藏Y轴刻度线
        rangeAxis.setAxisLineVisible(false);     // 隐藏Y轴线

        // 设置固定大小，禁用所有拖动和缩放
        plot.setDomainPannable(true);  // 启用domain拖动（用于右键拖动）
        plot.setRangePannable(true);   // 启用range拖动（用于右键拖动）

        // 禁止任何拉伸和缩放
        plot.setRangeZeroBaselineVisible(false);
        plot.setDomainZeroBaselineVisible(false);
        plot.setNoDataMessage("暂无数据");

        // 设置固定模式
        domainAxis.setAutoRangeStickyZero(false);
        rangeAxis.setAutoRangeStickyZero(false);

        // 设置固定坐标轴
        plot.setDomainAxis(domainAxis);
        plot.setRangeAxis(rangeAxis);

        // 设置渲染器
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setDefaultLinesVisible(false); // 只显示点，不显示线
        plot.setRenderer(renderer);

        // 设置标题字体
        chart.getTitle().setFont(new Font(Font.SANS_SERIF, Font.BOLD, 16));

        // 设置图例字体
        if (chart.getLegend() != null) {
            chart.getLegend().setItemFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        }
    }

    /**
     * 设置界面布局
     * 配置主窗口的组件布局和样式
     */
    private void setupLayout() {
        // 创建控制面板
        JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 10, 5));
        controlPanel.setBorder(BorderFactory.createTitledBorder("控制面板"));

        // 文件选择行
        JPanel filePanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        filePanel.add(new JLabel("文件:"));
        filePanel.add(filePathField);
        filePanel.add(selectFileButton);
        controlPanel.add(filePanel);

        // 参数输入行
        JPanel paramPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        paramPanel.add(new JLabel("eps:"));
        paramPanel.add(epsField);
        paramPanel.add(Box.createHorizontalStrut(10));
        paramPanel.add(new JLabel("minPts:"));
        paramPanel.add(minPtsField);
        paramPanel.add(Box.createHorizontalStrut(10));
        paramPanel.add(clusterButton);
        controlPanel.add(paramPanel);

        // 主布局
        setLayout(new BorderLayout(10, 10));
        add(controlPanel, BorderLayout.NORTH);

        // 自定义图例面板已在initializeComponents中初始化

        // 创建主内容面板，使用水平分割
        JPanel mainContentPanel = new JPanel(new BorderLayout());
        mainContentPanel.add(chartPanel, BorderLayout.CENTER);
        mainContentPanel.add(legendPanel, BorderLayout.EAST);

        add(mainContentPanel, BorderLayout.CENTER);

        // 添加边框
        ((JComponent) getContentPane()).setBorder(new EmptyBorder(10, 10, 10, 10));
    }

    /**
     * 配置统一鼠标处理器
     * 整合右键拖拽平移和左键框选放大功能，避免事件冲突
     *
     * @param chartPanel 图表面板组件
     */
    private void setupUnifiedMouseHandler(ChartPanel chartPanel) {
        // 右键拖拽相关变量
        final java.awt.Point[] dragStartPoint = new java.awt.Point[1]; // 拖拽起始点
        final double[] dragStartDomain = new double[2]; // X轴范围 [min, max]
        final double[] dragStartRange = new double[2];  // Y轴范围 [min, max]
        final boolean[] isDragging = new boolean[1]; // 是否正在拖拽

        // 左键框选相关变量
        final java.awt.Point[] selectionStartPoint = new java.awt.Point[1]; // 框选起始点
        final java.awt.Point[] selectionEndPoint = new java.awt.Point[1];   // 框选结束点
        final boolean[] isSelecting = new boolean[1]; // 是否正在框选
        final Rectangle[] selectionRectangle = new Rectangle[1]; // 框选矩形

        // 创建统一的鼠标处理器
        java.awt.event.MouseAdapter mouseAdapter = new java.awt.event.MouseAdapter() {
            @Override
            public void mousePressed(java.awt.event.MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e)) {
                    // 右键拖拽开始
                    dragStartPoint[0] = e.getPoint();
                    XYPlot plot = chart.getXYPlot();
                    dragStartDomain[0] = plot.getDomainAxis().getLowerBound();
                    dragStartDomain[1] = plot.getDomainAxis().getUpperBound();
                    dragStartRange[0] = plot.getRangeAxis().getLowerBound();
                    dragStartRange[1] = plot.getRangeAxis().getUpperBound();
                    isDragging[0] = true;
                    chartPanel.setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
                    e.consume(); // 消费事件，防止冒泡
                } else if (SwingUtilities.isLeftMouseButton(e) && coordinates.size() > 0) {
                    // 左键框选开始
                    selectionStartPoint[0] = e.getPoint();
                    selectionEndPoint[0] = e.getPoint();
                    isSelecting[0] = true;
                    chartPanel.setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
                    e.consume(); // 消费事件，防止冒泡
                }
            }

            @Override
            public void mouseReleased(java.awt.event.MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e) && isDragging[0]) {
                    // 右键拖拽结束
                    isDragging[0] = false;
                    chartPanel.setCursor(Cursor.getDefaultCursor());
                    e.consume(); // 消费事件，防止冒泡
                } else if (SwingUtilities.isLeftMouseButton(e) && isSelecting[0] && selectionStartPoint[0] != null) {
                    // 左键框选结束，执行放大
                    isSelecting[0] = false;
                    chartPanel.setCursor(Cursor.getDefaultCursor());

                    // 使用实际绘制的正方形框选区域坐标
                    int x = selectionRectangle[0].x;
                    int y = selectionRectangle[0].y;
                    int width = selectionRectangle[0].width;
                    int height = selectionRectangle[0].height;

                    // 确保框选区域足够大（避免误操作）
                    if (width > 10 && height > 10) {
                        // 将屏幕坐标转换为数据坐标
                        Rectangle2D dataArea = chartPanel.getScreenDataArea();
                        if (dataArea != null) {
                            XYPlot plot = chart.getXYPlot();
                            org.jfree.chart.axis.ValueAxis domainAxis = plot.getDomainAxis();
                            org.jfree.chart.axis.ValueAxis rangeAxis = plot.getRangeAxis();

                            // 计算框选区域的数据坐标范围
                            double domainMin = domainAxis.java2DToValue(x, dataArea, plot.getDomainAxisEdge());
                            double domainMax = domainAxis.java2DToValue(x + width, dataArea, plot.getDomainAxisEdge());
                            double rangeMin = rangeAxis.java2DToValue(y + height, dataArea, plot.getRangeAxisEdge()); // 注意Y轴方向
                            double rangeMax = rangeAxis.java2DToValue(y, dataArea, plot.getRangeAxisEdge()); // 注意Y轴方向

                            // 确保坐标顺序正确
                            if (domainMin > domainMax) {
                                double temp = domainMin;
                                domainMin = domainMax;
                                domainMax = temp;
                            }
                            if (rangeMin > rangeMax) {
                                double temp = rangeMin;
                                rangeMin = rangeMax;
                                rangeMax = temp;
                            }

                            // 计算框选区域和图表区域的比例，保持XY比例锁定
                            double selectedDomainRange = domainMax - domainMin;
                            double selectedRangeRange = rangeMax - rangeMin;
                            double dataAspectRatio = selectedDomainRange / selectedRangeRange;
                            double panelAspectRatio = dataArea.getWidth() / dataArea.getHeight();

                            // 调整范围以保持比例，模拟地图效果
                            if (dataAspectRatio > panelAspectRatio) {
                                // 域范围相对较大，需要扩展范围范围
                                double newRangeRange = selectedDomainRange / panelAspectRatio;
                                double rangeCenter = (rangeMin + rangeMax) / 2;
                                rangeMin = rangeCenter - newRangeRange / 2;
                                rangeMax = rangeCenter + newRangeRange / 2;
                            } else {
                                // 范围范围相对较大，需要扩展域范围
                                double newDomainRange = selectedRangeRange * panelAspectRatio;
                                double domainCenter = (domainMin + domainMax) / 2;
                                domainMin = domainCenter - newDomainRange / 2;
                                domainMax = domainCenter + newDomainRange / 2;
                            }

                            // 应用放大到框选区域，充满整个窗口
                            plot.getDomainAxis().setRange(domainMin, domainMax);
                            plot.getRangeAxis().setRange(rangeMin, rangeMax);
                        }
                    }

                    // 清除框选矩形
                    if (selectionRectangle[0] != null) {
                        Graphics g = chartPanel.getGraphics();
                        if (g != null) {
                            g.setXORMode(chartPanel.getBackground());
                            g.drawRect(selectionRectangle[0].x, selectionRectangle[0].y, selectionRectangle[0].width, selectionRectangle[0].height);
                            g.setPaintMode();
                        }
                        selectionRectangle[0] = null;
                    }
                    e.consume(); // 消费事件，防止冒泡
                }
            }
        };

        // 创建统一的鼠标移动处理器
        java.awt.event.MouseMotionAdapter mouseMotionAdapter = new java.awt.event.MouseMotionAdapter() {
            @Override
            public void mouseDragged(java.awt.event.MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e) && isDragging[0] && dragStartPoint[0] != null) {
                    // 右键拖拽处理
                    java.awt.Point currentPoint = e.getPoint();

                    // 计算像素移动距离
                    int deltaX = currentPoint.x - dragStartPoint[0].x;
                    int deltaY = currentPoint.y - dragStartPoint[0].y;

                    // 获取图表尺寸
                    int chartWidth = chartPanel.getWidth();
                    int chartHeight = chartPanel.getHeight();

                    if (chartWidth > 0 && chartHeight > 0) {
                        // 计算坐标轴范围
                        double domainRange = dragStartDomain[1] - dragStartDomain[0];
                        double rangeRange = dragStartRange[1] - dragStartRange[0];

                        // 将像素移动转换为坐标值移动
                        double domainDelta = -deltaX * (domainRange / chartWidth); // X轴坐标变化量
                        double rangeDelta = deltaY * (rangeRange / chartHeight); // Y轴坐标变化量

                        // 更新坐标轴范围（基于初始位置计算，避免累积误差）
                        XYPlot plot = chart.getXYPlot();
                        plot.getDomainAxis().setRange(dragStartDomain[0] + domainDelta, dragStartDomain[1] + domainDelta); // 更新X轴范围
                        plot.getRangeAxis().setRange(dragStartRange[0] + rangeDelta, dragStartRange[1] + rangeDelta); // 更新Y轴范围
                    }
                    e.consume(); // 消费事件，防止冒泡
                } else if (SwingUtilities.isLeftMouseButton(e) && isSelecting[0] && selectionStartPoint[0] != null) {
                    // 左键框选拖拽处理
                    selectionEndPoint[0] = e.getPoint();

                    // 清除之前的框选矩形
                    if (selectionRectangle[0] != null) {
                        Graphics g = chartPanel.getGraphics();
                        if (g != null) {
                            g.setXORMode(chartPanel.getBackground());
                            g.drawRect(selectionRectangle[0].x, selectionRectangle[0].y, selectionRectangle[0].width, selectionRectangle[0].height);
                            g.setPaintMode();
                        }
                    }

                    // 绘制新的框选矩形，强制保持正方形比例
                    int x1 = selectionStartPoint[0].x;
                    int y1 = selectionStartPoint[0].y;
                    int x2 = selectionEndPoint[0].x;
                    int y2 = selectionEndPoint[0].y;

                    // 计算鼠标移动的最大距离，作为正方形的边长
                    int distance = Math.min(Math.abs(x2 - x1), Math.abs(y2 - y1));

                    // 确定正方形的坐标（保持与鼠标移动方向一致）
                    int x = x1 < x2 ? x1 : x1 - distance;
                    int y = y1 < y2 ? y1 : y1 - distance;
                    int size = distance;

                    selectionRectangle[0] = new Rectangle(x, y, size, size);

                    // 绘制框选矩形
                    Graphics g = chartPanel.getGraphics();
                    if (g != null) {
                        g.setXORMode(Color.WHITE);
                        g.drawRect(x, y, size, size);
                        g.setPaintMode();
                    }
                    e.consume(); // 消费事件，防止冒泡
                }
            }
        };

        // 添加统一的鼠标处理器
        chartPanel.addMouseListener(mouseAdapter);
        chartPanel.addMouseMotionListener(mouseMotionAdapter);
    }

    /**
     * 配置鼠标滚轮缩放功能
     * 实现以鼠标位置为中心的滚轮缩放功能
     *
     * @param chartPanel 图表面板组件
     */
    private void setupMouseHoverZoom(ChartPanel chartPanel) {
        final java.awt.Point hoverPoint = new java.awt.Point();

        /**
         * 处理鼠标移动事件
         * 记录当前鼠标位置，用于后续滚轮缩放操作
         *
         * @param e 鼠标事件对象，包含当前鼠标位置信息
         */
        chartPanel.addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
            @Override
            public void mouseMoved(java.awt.event.MouseEvent e) {
                hoverPoint.setLocation(e.getPoint());
            }
        });
        /**
         * 处理鼠标滚轮事件
         * 实现以鼠标位置为中心的滚轮缩放功能
         *
         * @param e 鼠标滚轮事件对象，包含滚轮旋转信息
         */
        chartPanel.addMouseWheelListener(new java.awt.event.MouseWheelListener() {
            @Override
            public void mouseWheelMoved(java.awt.event.MouseWheelEvent e) {
                if (coordinates.isEmpty()) {
                    return;
                }

                // 将鼠标位置转换为数据坐标
                java.awt.geom.Point2D plotPoint = chartPanel.translateScreenToJava2D(hoverPoint);
                XYPlot plot = chart.getXYPlot();

                // 获取坐标轴
                org.jfree.chart.axis.ValueAxis domainAxis = plot.getDomainAxis();
                org.jfree.chart.axis.ValueAxis rangeAxis = plot.getRangeAxis();

                // 计算鼠标位置对应的数据坐标
                Rectangle2D dataArea = chartPanel.getScreenDataArea();
                if (dataArea == null) {
                    return;
                }

                double mouseX = domainAxis.java2DToValue(plotPoint.getX(), dataArea, plot.getDomainAxisEdge());
                double mouseY = rangeAxis.java2DToValue(plotPoint.getY(), dataArea, plot.getRangeAxisEdge());

                // 获取当前显示范围
                double currentDomainMin = domainAxis.getLowerBound();
                double currentDomainMax = domainAxis.getUpperBound();
                double currentRangeMin = rangeAxis.getLowerBound();
                double currentRangeMax = rangeAxis.getUpperBound();

                // 计算缩放比例（滚轮向上放大，向下缩小）
                double zoomFactor = e.getWheelRotation() < 0 ? 0.9 : 1.1; // 0.9放大（滚轮向上），1.1缩小（滚轮向下）

                // 获取图表面板的宽高比
                double panelAspectRatio = dataArea.getWidth() / dataArea.getHeight();

                // 计算新的显示范围（以鼠标位置为中心进行等比缩放）
                double newDomainRange = (currentDomainMax - currentDomainMin) * zoomFactor;
                double newRangeRange = (currentRangeMax - currentRangeMin) * zoomFactor;

                // 根据面板比例调整缩放，保持地图比例
                double adjustedRangeRange = newDomainRange / panelAspectRatio;

                double newDomainMin = mouseX - (mouseX - currentDomainMin) * zoomFactor;
                double newDomainMax = mouseX + (currentDomainMax - mouseX) * zoomFactor;
                double newRangeMin = mouseY - (mouseY - currentRangeMin) * (adjustedRangeRange / (currentRangeMax - currentRangeMin));
                double newRangeMax = mouseY + (currentRangeMax - mouseY) * (adjustedRangeRange / (currentRangeMax - currentRangeMin));

                // 应用缩放
                plot.getDomainAxis().setRange(newDomainMin, newDomainMax);
                plot.getRangeAxis().setRange(newRangeMin, newRangeMax);
            }
        });
    }


    /**
     * 配置事件监听
     * 设置按钮点击事件和文件选择事件的处理逻辑
     */
    private void setupEventListeners() {
        // 文件选择按钮
        selectFileButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                selectAndLoadFile();
            }
        });

        // 聚类按钮
        clusterButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                performClustering();
            }
        });
    }

    /**
     * 文件选择和加载
     * 打开文件选择对话框并加载选中的坐标文件
     */
    private void selectAndLoadFile() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("选择坐标文件");
        fileChooser.setFileFilter(new FileNameExtensionFilter("坐标文件 (*.txt, *.csv)", "txt", "csv"));

        // 设置默认目录
        File testDir = new File("testFiles");
        if (testDir.exists() && testDir.isDirectory()) {
            fileChooser.setCurrentDirectory(testDir);
        }

        int result = fileChooser.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            File selectedFile = fileChooser.getSelectedFile();
            loadCoordinateFile(selectedFile);
        }
    }

    /**
     * 加载坐标数据文件
     * 解析文件内容并提取坐标数据，支持CSV和TXT格式
     *
     * @param file 要加载的坐标文件
     */
    private void loadCoordinateFile(File file) {
        try {
            log.info("开始加载坐标文件: {}", file.getAbsolutePath());
            coordinates.clear();
            dataset.removeAllSeries();

            List<String> lines = FileUtil.readLines(file, "UTF-8");

            for (String line : lines) {
                line = line.trim();
                if (StrUtil.isEmpty(line)) continue;

                String[] parts = line.split(",");
                if (parts.length >= 2) {
                    try {
                        double x = Double.parseDouble(parts[0].trim());
                        double y = Double.parseDouble(parts[1].trim());
                        coordinates.add(new double[]{x, y});
                    } catch (NumberFormatException e) {
                        log.warn("解析坐标失败: {}", line);
                    }
                }
            }

            if (coordinates.isEmpty()) {
                JOptionPane.showMessageDialog(this, "文件中没有找到有效的坐标数据！", "错误", JOptionPane.ERROR_MESSAGE);
                return;
            }

            // 显示原始数据
            displayOriginalData();

            filePathField.setText(file.getName());
            clusterButton.setEnabled(true);


            log.info("成功加载 {} 个坐标点", coordinates.size());

        } catch (Exception e) {
            log.error("加载文件失败", e);
            JOptionPane.showMessageDialog(this, "文件加载失败: " + e.getMessage(), "错误", JOptionPane.ERROR_MESSAGE);
        }
    }

    /**
     * 显示原始数据点
     * 在图表上显示加载的原始坐标数据
     */
    private void displayOriginalData() {
        XYSeries series = new XYSeries("原始数据");

        for (double[] coord : coordinates) {
            series.add(coord[0], coord[1]);
        }

        dataset.removeAllSeries();
        dataset.addSeries(series);

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, false);
        renderer.setSeriesShape(0, new java.awt.geom.Ellipse2D.Double(-0.5, -0.5, 1, 1)); // 1x1点
        renderer.setSeriesPaint(0, Color.BLACK);
        renderer.setSeriesVisibleInLegend(0, false);

        chart.getXYPlot().setRenderer(renderer);
        centerChart();
    }

    /**
     * 执行DBSCAN聚类
     * 获取用户输入的参数并执行聚类算法
     */
    private void performClustering() {
        if (coordinates.isEmpty()) {
            JOptionPane.showMessageDialog(this, "请先加载数据文件！", "警告", JOptionPane.WARNING_MESSAGE);
            return;
        }

        try {
            double eps = Double.parseDouble(epsField.getText());
            int minPts = Integer.parseInt(minPtsField.getText());

            if (eps <= 0 || minPts < 2) {
                JOptionPane.showMessageDialog(this, "请输入有效的参数 (eps > 0, minPts ≥ 2)", "参数错误", JOptionPane.ERROR_MESSAGE);
                return;
            }

            log.info("开始执行DBSCAN聚类: eps={}, minPts={}", eps, minPts);

            double[][] data = new double[coordinates.size()][2];
            for (int i = 0; i < coordinates.size(); i++) {
                data[i] = coordinates.get(i);
            }

            int[] labels = performDBSCAN(data, eps, minPts);
            displayClusteringResults(labels, eps, minPts);

        } catch (NumberFormatException e) {
            JOptionPane.showMessageDialog(this, "参数格式错误，请输入有效的数字！", "参数错误", JOptionPane.ERROR_MESSAGE);
        } catch (Exception e) {
            log.error("聚类执行失败", e);
            JOptionPane.showMessageDialog(this, "聚类执行失败: " + e.getMessage(), "错误", JOptionPane.ERROR_MESSAGE);
        }
    }

    /**
     * 使用ELKI库执行DBSCAN聚类
     *
     * @param data   输入数据点数组
     * @param eps    邻域半径参数
     * @param minPts 最小点数参数
     *
     * @return 聚类标签数组，-1表示噪声点
     */
    private int[] performDBSCAN(double[][] data, double eps, int minPts) {
        // 创建数据库
        ArrayAdapterDatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);
        Database db = new StaticArrayDatabase(dbc, null);
        db.initialize();

        // 获取关系
        Relation<DoubleVector> relation = db.getRelation(TypeUtil.DOUBLE_VECTOR_FIELD);

        // 创建DBSCAN算法
        DBSCAN<DoubleVector> dbscan = new DBSCAN<>(EuclideanDistance.STATIC, eps, minPts);

        // 执行聚类
        Clustering<Model> result = dbscan.run(relation);

        // 获取聚类结果
        int[] labels = new int[data.length];
        Arrays.fill(labels, -1); // 初始化为噪声点，-1表示噪声

        int clusterIndex = 0;
        for (Cluster<Model> cluster : result.getAllClusters()) {
            if (!cluster.isNoise()) {
                // 设置簇标签
                DBIDRange ids = (DBIDRange) relation.getDBIDs();
                for (DBIDIter iter = cluster.getIDs().iter(); iter.valid(); iter.advance()) {
                    int index = ids.getOffset(iter);
                    labels[index] = clusterIndex;
                }
                clusterIndex++;
            }
        }

        return labels;
    }

    /**
     * 更新图例面板
     * 根据聚类结果创建自定义图例，显示各簇的颜色和点数统计
     *
     * @param clusters    聚类结果映射
     * @param noisePoints 噪声点列表
     */
    private void updateCustomLegend(Map<Integer, List<double[]>> clusters, List<double[]> noisePoints) {
        legendPanel.removeAll();

        legendPanel.add(Box.createVerticalStrut(10));

        // 添加聚类图例项
        int clusterIndex = 0;
        for (Map.Entry<Integer, List<double[]>> entry : clusters.entrySet()) {
            int clusterId = entry.getKey();
            List<double[]> points = entry.getValue();

            JPanel itemPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
            itemPanel.setBackground(Color.WHITE);
            itemPanel.setMaximumSize(new Dimension(180, 25));

            // 16x16颜色指示器
            JPanel colorIndicator = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    g.setColor(CLUSTER_COLORS[clusterId % CLUSTER_COLORS.length]);
                    Graphics2D g2d = (Graphics2D) g;
                    g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                    g2d.fill(new java.awt.geom.Ellipse2D.Double(0, 0, 16, 16));
                }
            };
            colorIndicator.setPreferredSize(new Dimension(12, 12));
            colorIndicator.setBackground(Color.WHITE);

            JLabel label = new JLabel(String.format("簇 %d (%d点)", clusterId, points.size()));
            label.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));

            itemPanel.add(colorIndicator);
            itemPanel.add(label);

            legendPanel.add(itemPanel);
            legendPanel.add(Box.createVerticalStrut(5));

            clusterIndex++;
        }

        // 添加噪声点图例项
        if (!noisePoints.isEmpty()) {
            JPanel noisePanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
            noisePanel.setBackground(Color.WHITE);
            noisePanel.setMaximumSize(new Dimension(180, 25));

            // 噪声点颜色指示器
            JPanel noiseIndicator = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    g.setColor(NOISE_COLOR);
                    Graphics2D g2d = (Graphics2D) g;
                    g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                    g2d.fill(new java.awt.geom.Ellipse2D.Double(0, 0, 16, 16));
                }
            };
            noiseIndicator.setPreferredSize(new Dimension(12, 12));
            noiseIndicator.setBackground(Color.WHITE);

            JLabel noiseLabel = new JLabel(String.format("噪声 (%d点)", noisePoints.size()));
            noiseLabel.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));

            noisePanel.add(noiseIndicator);
            noisePanel.add(noiseLabel);

            legendPanel.add(noisePanel);
        }

        legendPanel.add(Box.createVerticalGlue());
        legendPanel.revalidate();
        legendPanel.repaint();
    }

    /**
     * 显示聚类结果
     * 在图表上显示聚类结果，包括各簇点和噪声点
     *
     * @param labels 聚类标签数组
     * @param eps    使用的邻域半径参数
     * @param minPts 使用的最小点数参数
     */
    private void displayClusteringResults(int[] labels, double eps, int minPts) {
        dataset.removeAllSeries();

        Map<Integer, List<double[]>> clusters = new HashMap<>();
        List<double[]> noisePoints = new ArrayList<>();

        for (int i = 0; i < labels.length; i++) {
            int label = labels[i];
            if (label == -1) {
                noisePoints.add(coordinates.get(i)); // 噪声点
            } else {
                clusters.computeIfAbsent(label, k -> new ArrayList<>()).add(coordinates.get(i)); // 聚类点
            }
        }

        int clusterCount = clusters.size();
        int noiseCount = noisePoints.size();

        updateCustomLegend(clusters, noisePoints);

        int seriesIndex = 0;

        // 添加聚类点
        for (Map.Entry<Integer, List<double[]>> entry : clusters.entrySet()) {
            int clusterId = entry.getKey();
            List<double[]> points = entry.getValue();

            XYSeries series = new XYSeries(String.format("簇 %d (%d点)", clusterId, points.size()));

            for (double[] point : points) {
                series.add(point[0], point[1]); // 添加坐标点到序列
            }

            dataset.addSeries(series);

            XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) chart.getXYPlot().getRenderer();
            renderer.setSeriesLinesVisible(seriesIndex, false); // 只显示点，不显示线
            renderer.setSeriesShape(seriesIndex, new java.awt.geom.Ellipse2D.Double(-0.5, -0.5, 1, 1)); // 1x1点
            renderer.setSeriesPaint(seriesIndex, CLUSTER_COLORS[clusterId % CLUSTER_COLORS.length]); // 设置聚类颜色

            seriesIndex++;
        }

        // 添加噪声点
        if (!noisePoints.isEmpty()) {
            XYSeries noiseSeries = new XYSeries(String.format("噪声 (%d点)", noiseCount));

            for (double[] point : noisePoints) {
                noiseSeries.add(point[0], point[1]); // 添加噪声点到序列
            }

            dataset.addSeries(noiseSeries);

            XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) chart.getXYPlot().getRenderer();
            renderer.setSeriesLinesVisible(seriesIndex, false); // 只显示点，不显示线
            renderer.setSeriesShape(seriesIndex, new java.awt.geom.Ellipse2D.Double(-0.5, -0.5, 1, 1)); // 1x1点
            renderer.setSeriesPaint(seriesIndex, NOISE_COLOR); // 设置噪声点颜色
        }

        centerChart(); // 居中显示所有数据点

        log.info("聚类完成: {}个簇, {}个噪声点, 共{}个点", clusterCount, noiseCount, coordinates.size()); // 记录聚类结果统计
    }


    /**
     * 获取X坐标最小值
     *
     * @return X坐标最小值
     */
    private double getMinX() {
        return coordinates.stream().mapToDouble(c -> c[0]).min().orElse(0);
    }

    /**
     * 获取X坐标最大值
     *
     * @return X坐标最大值
     */
    private double getMaxX() {
        return coordinates.stream().mapToDouble(c -> c[0]).max().orElse(0);
    }

    /**
     * 获取Y坐标最小值
     *
     * @return Y坐标最小值
     */
    private double getMinY() {
        return coordinates.stream().mapToDouble(c -> c[1]).min().orElse(0);
    }

    /**
     * 获取Y坐标最大值
     *
     * @return Y坐标最大值
     */
    private double getMaxY() {
        return coordinates.stream().mapToDouble(c -> c[1]).max().orElse(0);
    }

    /**
     * 居中显示所有数据点
     * 自动调整图表显示范围，确保所有数据点都能在视图中完整显示
     */
    private void centerChart() {
        if (coordinates.isEmpty()) {
            return;
        }

        double minX = getMinX();
        double maxX = getMaxX();
        double minY = getMinY();
        double maxY = getMaxY();

        // 计算数据范围
        double xRange = maxX - minX;
        double yRange = maxY - minY;

        // 获取图表面板的宽高比
        int panelWidth = chartPanel.getWidth();
        int panelHeight = chartPanel.getHeight();

        if (panelWidth <= 0 || panelHeight <= 0) {
            panelWidth = 800;  // 默认值
            panelHeight = 600; // 默认值
        }

        // 获取图表面板的宽高比
        double panelAspectRatio = (double) panelWidth / panelHeight; // 面板宽高比

        double dataAspectRatio = xRange / yRange; // 数据宽高比
        double margin = Math.max(xRange, yRange) * 0.1; // 边距为数据范围的10%

        double xMin, xMax, yMin, yMax;

        if (dataAspectRatio > panelAspectRatio) {
            // 数据比面板更宽，以X轴为主进行适配
            xMin = minX - margin;
            xMax = maxX + margin;
            double adjustedYRange = (xMax - xMin) / panelAspectRatio; // 根据面板比例调整Y轴范围
            double yCenter = (minY + maxY) / 2; // Y轴中心点
            yMin = yCenter - adjustedYRange / 2;
            yMax = yCenter + adjustedYRange / 2;
        } else {
            // 数据比面板更高，以Y轴为主进行适配
            yMin = minY - margin;
            yMax = maxY + margin;
            double adjustedXRange = (yMax - yMin) * panelAspectRatio; // 根据面板比例调整X轴范围
            double xCenter = (minX + maxX) / 2; // X轴中心点
            xMin = xCenter - adjustedXRange / 2;
            xMax = xCenter + adjustedXRange / 2;
        }

        // 设置坐标轴范围
        XYPlot plot = chart.getXYPlot();
        plot.getDomainAxis().setRange(xMin, xMax); // 设置X轴范围
        plot.getRangeAxis().setRange(yMin, yMax); // 设置Y轴范围
    }

    /**
     * 程序入口
     * 启动DBSCAN聚类可视化应用程序
     *
     * @param args 命令行参数
     */
    public static void main(String[] args) {
        try {
            // 设置系统外观为默认外观，确保在不同操作系统上有一致的外观
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            log.warn("设置系统外观失败", e);
        }

        /**
         * 启动Swing事件调度线程
         * 确保GUI组件在事件线程中创建和更新
         */
        SwingUtilities.invokeLater(() -> {
            ClusterVisualizationGui gui = new ClusterVisualizationGui();
            gui.setVisible(true);

            // 窗口大小变化时重新居中
            gui.addComponentListener(new ComponentAdapter() {
                @Override
                public void componentResized(ComponentEvent e) {
                    if (gui.coordinates != null && !gui.coordinates.isEmpty()) {
                        gui.centerChart();
                    }
                }
            });
        });
    }
}