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
import org.jfree.chart.ui.HorizontalAlignment;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.VerticalAlignment;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
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
 * DBSCAN聚类可视化GUI工具
 * 模仿Python版本的功能，使用Swing + JFreeChart实现
 *
 * @author sunyu
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

    // 颜色配置 - 使用更鲜艳的颜色
    private static final Color[] CLUSTER_COLORS = {new Color(255, 99, 71),   // 番茄红
            new Color(30, 144, 255),  // 道奇蓝
            new Color(50, 205, 50),   // 酸橙绿
            new Color(255, 215, 0),   // 金色
            new Color(138, 43, 226),  // 蓝紫色
            new Color(255, 140, 0),   // 深橙色
            new Color(220, 20, 60),   // 深粉红色
            new Color(0, 191, 255),   // 深天蓝
            new Color(255, 20, 147),  // 深粉红色
            new Color(0, 250, 154),   // 中等春绿色
            new Color(255, 69, 0),    // 橙红色
            new Color(106, 90, 205),  // 石板蓝
            new Color(255, 255, 0),   // 黄色
            new Color(186, 85, 211),  // 中等紫色
            new Color(0, 255, 127),   // 春绿色
            new Color(123, 104, 238)  // 中等石板蓝
    };
    private static final Color NOISE_COLOR = Color.BLACK;

    public ClusterVisualizationGui() {
        initializeComponents();
        setupLayout();
        setupEventListeners();
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setTitle("DBSCAN聚类可视化工具");
        setSize(1200, 800);
        setLocationRelativeTo(null); // 居中显示
    }

    /**
     * 初始化GUI组件
     */
    private void initializeComponents() {
        // 文件选择组件
        filePathField = new JTextField(100);
        filePathField.setEditable(false);
        filePathField.setToolTipText("选择的坐标文件路径");

        selectFileButton = new JButton("选择数据文件");
        selectFileButton.setToolTipText("选择TXT或CSV格式的坐标文件");

        // 聚类参数组件
        epsField = new JTextField("5", 8);
        epsField.setToolTipText("DBSCAN eps参数 (邻域半径)");

        minPtsField = new JTextField("20", 8);
        minPtsField.setToolTipText("DBSCAN minPts参数 (最小点数)");

        clusterButton = new JButton("执行聚类");
        clusterButton.setEnabled(false);
        clusterButton.setToolTipText("执行DBSCAN聚类分析");

        // 图表组件
        dataset = new XYSeriesCollection();
        chart = ChartFactory.createScatterPlot("DBSCAN聚类分析", "", "", dataset, PlotOrientation.VERTICAL, true,  // 显示图例
                true,  // 工具提示
                false  // URL链接
        );

        // 自定义图表样式
        customizeChart();
        
        // 隐藏默认图例
        chart.removeLegend();
        
        // 初始化自定义图例面板
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

        chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        chartPanel.setMouseWheelEnabled(false); // 禁用滚轮缩放
        chartPanel.setMouseZoomable(false);     // 禁用鼠标缩放
        chartPanel.setFillZoomRectangle(false);
        // 设置固定大小，不允许任何缩放
        chartPanel.setMaximumDrawWidth(800);
        chartPanel.setMaximumDrawHeight(600);
        chartPanel.setMinimumDrawWidth(800);
        chartPanel.setMinimumDrawHeight(600);

        // 完全禁用缩放和拖动功能
        chartPanel.setDomainZoomable(false);   // 禁用domain轴缩放
        chartPanel.setRangeZoomable(false);    // 禁用range轴缩放
        chartPanel.setPopupMenu(null);         // 禁用右键菜单
        
        // 设置等比例显示，保持地图比例
        chartPanel.setRefreshBuffer(true);
        chartPanel.setDoubleBuffered(true);

        // 添加右键拖动功能
        setupRightClickDrag(chartPanel);
        
        // 设置鼠标滚轮放大缩小功能
        setupMouseHoverZoom(chartPanel);
    }

    /**
     * 自定义图表样式
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
        plot.setDomainGridlinesVisible(true);
        plot.setRangeGridlinesVisible(true);

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
     * 设置布局
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
     * 设置右键拖动功能
     */
    private void setupRightClickDrag(ChartPanel chartPanel) {
        // 使用共享变量来避免重复初始化
        final java.awt.Point[] dragStartPoint = new java.awt.Point[1];
        final double[] dragStartDomain = new double[2]; // [min, max]
        final double[] dragStartRange = new double[2];  // [min, max]
        final boolean[] isDragging = new boolean[1];

        chartPanel.addMouseListener(new java.awt.event.MouseAdapter() {
            @Override
            public void mousePressed(java.awt.event.MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e)) {
                    dragStartPoint[0] = e.getPoint();
                    XYPlot plot = chart.getXYPlot();
                    dragStartDomain[0] = plot.getDomainAxis().getLowerBound();
                    dragStartDomain[1] = plot.getDomainAxis().getUpperBound();
                    dragStartRange[0] = plot.getRangeAxis().getLowerBound();
                    dragStartRange[1] = plot.getRangeAxis().getUpperBound();
                    isDragging[0] = true;
                    chartPanel.setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
                }
            }

            @Override
            public void mouseReleased(java.awt.event.MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e)) {
                    isDragging[0] = false;
                    chartPanel.setCursor(Cursor.getDefaultCursor());
                }
            }
        });

        chartPanel.addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
            @Override
            public void mouseDragged(java.awt.event.MouseEvent e) {
                if (SwingUtilities.isRightMouseButton(e) && isDragging[0] && dragStartPoint[0] != null) {
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
                        double domainDelta = -deltaX * (domainRange / chartWidth);
                        double rangeDelta = deltaY * (rangeRange / chartHeight);

                        // 更新坐标轴范围（基于初始位置计算，避免累积误差）
                        XYPlot plot = chart.getXYPlot();
                        plot.getDomainAxis().setRange(dragStartDomain[0] + domainDelta, dragStartDomain[1] + domainDelta);
                        plot.getRangeAxis().setRange(dragStartRange[0] + rangeDelta, dragStartRange[1] + rangeDelta);
                    }
                }
            }
        });
    }

    /**
     * 设置鼠标滚轮放大缩小功能，保持地图比例
     */
    private void setupMouseHoverZoom(ChartPanel chartPanel) {
        final java.awt.Point hoverPoint = new java.awt.Point();
        
        chartPanel.addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
            @Override
            public void mouseMoved(java.awt.event.MouseEvent e) {
                hoverPoint.setLocation(e.getPoint());
            }
        });
        
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
                double zoomFactor = e.getWheelRotation() < 0 ? 0.9 : 1.1; // 0.9放大，1.1缩小
                
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
                
                // 限制最大缩放范围（防止缩放到无限大）
                double dataRange = Math.max(getMaxX() - getMinX(), getMaxY() - getMinY());
                double maxRange = dataRange * 100; // 最大100倍数据范围
                double minRange = dataRange * 0.01; // 最小1%数据范围
                
                if (newDomainRange >= minRange && newDomainRange <= maxRange && 
                    newRangeRange >= minRange && newRangeRange <= maxRange) {
                    // 应用缩放
                    plot.getDomainAxis().setRange(newDomainMin, newDomainMax);
                    plot.getRangeAxis().setRange(newRangeMin, newRangeMax);
                }
                
                // 应用缩放
                plot.getDomainAxis().setRange(newDomainMin, newDomainMax);
                plot.getRangeAxis().setRange(newRangeMin, newRangeMax);
            }
        });
    }

    /**
     * 设置事件监听器
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
     * 选择并加载文件
     */
    private void selectAndLoadFile() {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("选择坐标文件");

        // 设置文件过滤器
        FileNameExtensionFilter filter = new FileNameExtensionFilter("坐标文件 (*.txt, *.csv)", "txt", "csv");
        fileChooser.setFileFilter(filter);

        // 设置默认目录（如果有测试数据）
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
     * 加载坐标文件
     */
    private void loadCoordinateFile(File file) {
        try {
            log.info("开始加载坐标文件: {}", file.getAbsolutePath());
            coordinates.clear();
            dataset.removeAllSeries();

            // 读取文件内容
            List<String> lines = FileUtil.readLines(file, "UTF-8");

            for (String line : lines) {
                line = line.trim();
                if (StrUtil.isEmpty(line)) {
                    continue;
                }

                // 解析坐标
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
     * 显示原始数据
     */
    private void displayOriginalData() {
        XYSeries series = new XYSeries("原始数据");

        for (double[] coord : coordinates) {
            series.add(coord[0], coord[1]);
        }

        dataset.removeAllSeries();
        dataset.addSeries(series);

        // 设置点的样式
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, false);
        renderer.setSeriesShape(0, new java.awt.geom.Ellipse2D.Double(-1, -1, 2, 2));
        renderer.setSeriesPaint(0, Color.BLACK);
        renderer.setSeriesVisibleInLegend(0, false);

        chart.getXYPlot().setRenderer(renderer);
        chart.setTitle("原始数据 - 请选择文件并点击'执行聚类'按钮");

        // 居中显示所有数据点
        centerChart();
    }

    /**
     * 执行聚类
     */
    private void performClustering() {
        if (coordinates.isEmpty()) {
            JOptionPane.showMessageDialog(this, "请先加载数据文件！", "警告", JOptionPane.WARNING_MESSAGE);
            return;
        }

        try {
            // 获取参数
            double eps = Double.parseDouble(epsField.getText());
            int minPts = Integer.parseInt(minPtsField.getText());

            // 参数验证
            if (eps <= 0 || minPts < 2) {
                JOptionPane.showMessageDialog(this, "请输入有效的参数 (eps > 0, minPts ≥ 2)", "参数错误", JOptionPane.ERROR_MESSAGE);
                return;
            }

            log.info("开始执行DBSCAN聚类: eps={}, minPts={}", eps, minPts);

            // 准备数据用于ELKI
            double[][] data = new double[coordinates.size()][2];
            for (int i = 0; i < coordinates.size(); i++) {
                data[i] = coordinates.get(i);
            }

            // 执行DBSCAN聚类
            int[] labels = performDBSCAN(data, eps, minPts);

            // 显示聚类结果
        displayClusteringResults(labels, eps, minPts);

    } catch (NumberFormatException e) {
        JOptionPane.showMessageDialog(this, "参数格式错误，请输入有效的数字！", "参数错误", JOptionPane.ERROR_MESSAGE);
    } catch (Exception e) {
        log.error("聚类执行失败", e);
        JOptionPane.showMessageDialog(this, "聚类执行失败: " + e.getMessage(), "错误", JOptionPane.ERROR_MESSAGE);
    }
}

    /**
     * 使用ELKI执行DBSCAN聚类
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
        Arrays.fill(labels, -1); // 初始化为噪声点

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
     * 更新自定义图例面板
     */
    private void updateCustomLegend(Map<Integer, List<double[]>> clusters, List<double[]> noisePoints) {
        legendPanel.removeAll();
        
        // 添加标题
        JLabel titleLabel = new JLabel("聚类图例");
        titleLabel.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 14));
        titleLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        legendPanel.add(titleLabel);
        legendPanel.add(Box.createVerticalStrut(10));
        
        // 添加聚类图例项
        int clusterIndex = 0;
        for (Map.Entry<Integer, List<double[]>> entry : clusters.entrySet()) {
            int clusterId = entry.getKey();
            List<double[]> points = entry.getValue();
            
            JPanel itemPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
            itemPanel.setBackground(Color.WHITE);
            itemPanel.setMaximumSize(new Dimension(180, 25));
            
            // 创建颜色指示器（4x4椭圆）
            JPanel colorIndicator = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    g.setColor(CLUSTER_COLORS[clusterId % CLUSTER_COLORS.length]);
                    Graphics2D g2d = (Graphics2D) g;
                    g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                    g2d.fill(new java.awt.geom.Ellipse2D.Double(2, 2, 8, 8));
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
            
            // 创建噪声点颜色指示器（4x4椭圆，黑色）
            JPanel noiseIndicator = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    g.setColor(NOISE_COLOR);
                    Graphics2D g2d = (Graphics2D) g;
                    g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                    g2d.fill(new java.awt.geom.Ellipse2D.Double(2, 2, 8, 8));
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
        
        // 添加弹性空间
        legendPanel.add(Box.createVerticalGlue());
        
        legendPanel.revalidate();
        legendPanel.repaint();
    }

    /**
     * 显示聚类结果
     */
    private void displayClusteringResults(int[] labels, double eps, int minPts) {
        dataset.removeAllSeries();

        // 统计聚类信息
        Map<Integer, List<double[]>> clusters = new HashMap<>();
        List<double[]> noisePoints = new ArrayList<>();

        for (int i = 0; i < labels.length; i++) {
            int label = labels[i];
            if (label == -1) {
                noisePoints.add(coordinates.get(i));
            } else {
                clusters.computeIfAbsent(label, k -> new ArrayList<>()).add(coordinates.get(i));
            }
        }

        int clusterCount = clusters.size();
        int noiseCount = noisePoints.size();
        
        // 更新自定义图例
        updateCustomLegend(clusters, noisePoints);

        // 创建系列数据
        int seriesIndex = 0;

        // 添加聚类点
        for (Map.Entry<Integer, List<double[]>> entry : clusters.entrySet()) {
            int clusterId = entry.getKey();
            List<double[]> points = entry.getValue();

            XYSeries series = new XYSeries(String.format("簇 %d (%d点)", clusterId, points.size()));

            for (double[] point : points) {
                series.add(point[0], point[1]);
            }

            dataset.addSeries(series);

            // 设置系列样式
            XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) chart.getXYPlot().getRenderer();
            renderer.setSeriesLinesVisible(seriesIndex, false);
            renderer.setSeriesShape(seriesIndex, new java.awt.geom.Ellipse2D.Double(-0.5, -0.5, 1, 1));
            renderer.setSeriesPaint(seriesIndex, CLUSTER_COLORS[clusterId % CLUSTER_COLORS.length]);

            seriesIndex++;
        }

        // 添加噪声点
        if (!noisePoints.isEmpty()) {
            XYSeries noiseSeries = new XYSeries(String.format("噪声 (%d点)", noiseCount));

            for (double[] point : noisePoints) {
                noiseSeries.add(point[0], point[1]);
            }

            dataset.addSeries(noiseSeries);

            XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) chart.getXYPlot().getRenderer();
            renderer.setSeriesLinesVisible(seriesIndex, false);
            renderer.setSeriesShape(seriesIndex, new java.awt.geom.Ellipse2D.Double(-1.5, -1.5, 3, 3));
            renderer.setSeriesPaint(seriesIndex, NOISE_COLOR);
        }

        // 更新标题和统计信息
        chart.setTitle("DBSCAN聚类结果 (eps=" + eps + ", minPts=" + minPts + ")");

        String stats = "聚类统计\n" + "簇数量: " + clusterCount + "\n" + "噪声点: " + noiseCount + " (" + String.format("%.1f", (double) noiseCount / coordinates.size() * 100) + "%)\n" + "总点数: " + coordinates.size() + "\n" + "参数: eps=" + eps + ", minPts=" + minPts;

        // 居中显示所有数据点
        centerChart();

        log.info("聚类完成: {}个簇, {}个噪声点, 共{}个点", clusterCount, noiseCount, coordinates.size());
    }



    /**
     * 获取坐标范围
     */
    private double getMinX() {
        return coordinates.stream().mapToDouble(c -> c[0]).min().orElse(0);
    }

    private double getMaxX() {
        return coordinates.stream().mapToDouble(c -> c[0]).max().orElse(0);
    }

    private double getMinY() {
        return coordinates.stream().mapToDouble(c -> c[1]).min().orElse(0);
    }

    private double getMaxY() {
        return coordinates.stream().mapToDouble(c -> c[1]).max().orElse(0);
    }

    /**
     * 将图表居中显示所有数据点，保持地图比例
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
        
        double panelAspectRatio = (double) panelWidth / panelHeight;
        
        // 计算数据的高斯投影比例（通常接近1:1，但需要根据实际数据调整）
        double dataAspectRatio = xRange / yRange;
        
        // 添加边距
        double margin = Math.max(xRange, yRange) * 0.1;
        
        double xMin, xMax, yMin, yMax;
        
        if (dataAspectRatio > panelAspectRatio) {
            // 数据更宽，以X轴为主，调整Y轴范围
            xMin = minX - margin;
            xMax = maxX + margin;
            
            double adjustedYRange = (xMax - xMin) / panelAspectRatio;
            double yCenter = (minY + maxY) / 2;
            yMin = yCenter - adjustedYRange / 2;
            yMax = yCenter + adjustedYRange / 2;
        } else {
            // 数据更高，以Y轴为主，调整X轴范围
            yMin = minY - margin;
            yMax = maxY + margin;
            
            double adjustedXRange = (yMax - yMin) * panelAspectRatio;
            double xCenter = (minX + maxX) / 2;
            xMin = xCenter - adjustedXRange / 2;
            xMax = xCenter + adjustedXRange / 2;
        }

        XYPlot plot = chart.getXYPlot();
        plot.getDomainAxis().setRange(xMin, xMax);
        plot.getRangeAxis().setRange(yMin, yMax);
    }

    /**
     * 主方法 - 启动应用
     */
    public static void main(String[] args) {
        // 设置系统外观
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            log.warn("设置系统外观失败", e);
        }

        // 在事件调度线程中启动GUI
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                ClusterVisualizationGui gui = new ClusterVisualizationGui();
                gui.setVisible(true);

                // 添加组件监听器，确保窗口大小变化时保持图表比例
                gui.addComponentListener(new ComponentAdapter() {
                    @Override
                    public void componentResized(ComponentEvent e) {
                        // 重新设置图表的等比显示
                        if (gui.coordinates != null && !gui.coordinates.isEmpty()) {
                            gui.centerChart();
                        }
                    }
                });
            }
        });
    }
}