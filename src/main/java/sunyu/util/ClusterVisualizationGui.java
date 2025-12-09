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
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
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
    private JTextArea statsTextArea;
    private JFreeChart chart;

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
        filePathField = new JTextField(50);
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
        chart = ChartFactory.createScatterPlot("DBSCAN聚类分析", "X 坐标 (米)", "Y 坐标 (米)", dataset, PlotOrientation.VERTICAL, true,  // 显示图例
                true,  // 工具提示
                false  // URL链接
        );

        // 自定义图表样式
        customizeChart();

        chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        chartPanel.setMouseWheelEnabled(true); // 启用滚轮缩放
        chartPanel.setMouseZoomable(true);     // 启用鼠标缩放

        // 统计信息组件
        statsTextArea = new JTextArea(6, 30);
        statsTextArea.setEditable(false);
        statsTextArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        statsTextArea.setBorder(BorderFactory.createTitledBorder("聚类统计信息"));
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

        // 设置坐标轴
        NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();

        domainAxis.setAutoRangeIncludesZero(false);
        rangeAxis.setAutoRangeIncludesZero(false);
        domainAxis.setLabelFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        rangeAxis.setLabelFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));
        
        // 设置坐标轴刻度字体
        domainAxis.setTickLabelFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));
        rangeAxis.setTickLabelFont(new Font(Font.SANS_SERIF, Font.PLAIN, 10));

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

        // 统计信息面板
        JScrollPane statsScrollPane = new JScrollPane(statsTextArea);
        statsScrollPane.setPreferredSize(new Dimension(300, 120));
        statsScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        // 主布局
        setLayout(new BorderLayout(10, 10));
        add(controlPanel, BorderLayout.NORTH);

        // 中间区域：图表和统计信息
        JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        splitPane.setLeftComponent(chartPanel);
        splitPane.setRightComponent(statsScrollPane);
        splitPane.setDividerLocation(900);
        splitPane.setOneTouchExpandable(true);

        add(splitPane, BorderLayout.CENTER);

        // 添加边框
        ((JComponent) getContentPane()).setBorder(new EmptyBorder(10, 10, 10, 10));
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

            // 更新统计信息
            updateStats("数据加载成功！\n" + "数据点数量: " + coordinates.size() + "\n" + "X坐标范围: [" + String.format("%.2f", getMinX()) + ", " + String.format("%.2f", getMaxX()) + "]\n" + "Y坐标范围: [" + String.format("%.2f", getMinY()) + ", " + String.format("%.2f", getMaxY()) + "]");

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
        renderer.setSeriesShape(0, new java.awt.geom.Ellipse2D.Double(-4, -4, 8, 8));
        renderer.setSeriesPaint(0, Color.BLACK);
        renderer.setSeriesVisibleInLegend(0, false);

        chart.getXYPlot().setRenderer(renderer);
        chart.setTitle("原始数据 - 请选择文件并点击'执行聚类'按钮");
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
            renderer.setSeriesShape(seriesIndex, new java.awt.geom.Ellipse2D.Double(-5, -5, 10, 10));
            renderer.setSeriesPaint(seriesIndex, CLUSTER_COLORS[clusterId % CLUSTER_COLORS.length]);
            renderer.setSeriesVisibleInLegend(seriesIndex, true);  // 确保簇系列在图例中显示

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
            renderer.setSeriesShape(seriesIndex, new java.awt.geom.Ellipse2D.Double(-4, -4, 8, 8));
            renderer.setSeriesPaint(seriesIndex, NOISE_COLOR);
            renderer.setSeriesVisibleInLegend(seriesIndex, true);  // 确保噪声点在图例中显示
        }

        // 更新标题和统计信息
        chart.setTitle("DBSCAN聚类结果 (eps=" + eps + ", minPts=" + minPts + ")");

        String stats = "聚类统计\n" + "簇数量: " + clusterCount + "\n" + "噪声点: " + noiseCount + " (" + String.format("%.1f", (double) noiseCount / coordinates.size() * 100) + "%)\n" + "总点数: " + coordinates.size() + "\n" + "参数: eps=" + eps + ", minPts=" + minPts;

        updateStats(stats);

        log.info("聚类完成: {}个簇, {}个噪声点, 共{}个点", clusterCount, noiseCount, coordinates.size());
    }

    /**
     * 更新统计信息
     */
    private void updateStats(String stats) {
        statsTextArea.setText(stats);
        statsTextArea.setCaretPosition(0); // 滚动到顶部
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

                // 显示欢迎信息
                gui.updateStats("操作步骤：\n" + "1. 点击'选择数据文件'加载坐标数据\n" + "2. 设置聚类参数 (eps和minPts)\n" + "3. 点击'执行聚类'进行分析\n\n" + "支持TXT和CSV格式的坐标文件，\n" + "每行格式：x,y");
            }
        });
    }
}