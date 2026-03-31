package sunyu.util.test;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;

import javax.swing.*;
import java.awt.*;
import java.awt.Point;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;

/**
 * WKT多边形可视化工具
 * 支持输入Polygon或MultiPolygon的WKT字符串，并在画布上直接绘制轮廓
 *
 * @author SunYu
 */
public class ShowGaussPointsGui extends JFrame {

    private JTextArea wktTextArea;
    private JButton displayButton;
    private JLabel statusLabel;
    private PolygonCanvas canvas;

    private final GeometryFactory geometryFactory = new GeometryFactory();
    private final WKTReader wktReader = new WKTReader(geometryFactory);

    public ShowGaussPointsGui() {
        initComponents();
        setupLayout();
        setupListeners();

        setTitle("WKT多边形可视化工具");
        setSize(1200, 800);
        setLocationRelativeTo(null);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    private void initComponents() {
        wktTextArea = new JTextArea(6, 80);
        wktTextArea.setLineWrap(true);
        wktTextArea.setWrapStyleWord(true);
        wktTextArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));

        displayButton = new JButton("显示轮廓");
        displayButton.setFont(new Font(Font.SANS_SERIF, Font.BOLD, 14));

        statusLabel = new JLabel("请输入WKT字符串，支持Polygon或MultiPolygon");
        statusLabel.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 12));

        canvas = new PolygonCanvas();
    }

    private void setupLayout() {
        setLayout(new BorderLayout(10, 10));

        JPanel topPanel = new JPanel(new BorderLayout(5, 5));
        topPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 5, 10));

        JPanel statusPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        statusPanel.add(statusLabel);
        topPanel.add(statusPanel, BorderLayout.NORTH);

        JScrollPane scrollPane = new JScrollPane(wktTextArea);
        scrollPane.setBorder(BorderFactory.createTitledBorder("WKT输入区域 (Polygon / MultiPolygon)"));
        topPanel.add(scrollPane, BorderLayout.CENTER);

        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        buttonPanel.add(displayButton);
        topPanel.add(buttonPanel, BorderLayout.SOUTH);

        add(topPanel, BorderLayout.NORTH);

        JPanel canvasContainer = new JPanel(new BorderLayout());
        canvasContainer.setBorder(BorderFactory.createEmptyBorder(5, 10, 10, 10));
        canvasContainer.add(canvas, BorderLayout.CENTER);
        add(canvasContainer, BorderLayout.CENTER);
    }

    private void setupListeners() {
        displayButton.addActionListener(e -> displayWkt());

        wktTextArea.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ENTER && e.isControlDown()) {
                    displayWkt();
                }
            }
        });
    }

    private void displayWkt() {
        String wktString = wktTextArea.getText().trim();

        if (wktString.isEmpty()) {
            statusLabel.setText("错误：WKT字符串不能为空");
            statusLabel.setForeground(Color.RED);
            return;
        }

        try {
            Geometry geometry = wktReader.read(wktString);

            if (!(geometry instanceof Polygon) && !(geometry instanceof MultiPolygon)) {
                statusLabel.setText("错误：WKT必须是Polygon或MultiPolygon类型，当前类型：" + geometry.getGeometryType());
                statusLabel.setForeground(Color.RED);
                return;
            }

            canvas.setGeometry(geometry);

            statusLabel.setText("成功显示：" + geometry.getGeometryType() + "，共" + geometry.getNumGeometries() + "个几何图形");
            statusLabel.setForeground(new Color(0, 128, 0));

        } catch (ParseException e) {
            statusLabel.setText("错误：WKT解析失败 - " + e.getMessage());
            statusLabel.setForeground(Color.RED);
        } catch (Exception e) {
            statusLabel.setText("错误：" + e.getMessage());
            statusLabel.setForeground(Color.RED);
        }
    }

    /**
     * 自定义画布面板，用于绘制多边形
     */
    private class PolygonCanvas extends JPanel {
        private Geometry geometry;
        private double viewCenterX = 0;  // 视图中心在几何坐标系中的X位置
        private double viewCenterY = 0;  // 视图中心在几何坐标系中的Y位置
        private double scale = 1.0;      // 缩放比例
        private Point lastDragPoint;
        private double autoScale = 1.0;  // 自动缩放比例，用于适应画布
        private double geomCenterX = 0;  // 几何图形中心X
        private double geomCenterY = 0;  // 几何图形中心Y

        public PolygonCanvas() {
            setBackground(Color.WHITE);
            setBorder(BorderFactory.createLineBorder(Color.GRAY));

            MouseAdapter mouseAdapter = new MouseAdapter() {
                @Override
                public void mousePressed(MouseEvent e) {
                    lastDragPoint = e.getPoint();
                }

                @Override
                public void mouseDragged(MouseEvent e) {
                    if (geometry == null) return;
                    Point currentPoint = e.getPoint();
                    // 将屏幕位移转换为几何坐标位移
                    double dx = (currentPoint.x - lastDragPoint.x) / (autoScale * scale);
                    double dy = -(currentPoint.y - lastDragPoint.y) / (autoScale * scale); // Y轴翻转
                    viewCenterX -= dx;
                    viewCenterY -= dy;
                    lastDragPoint = currentPoint;
                    repaint();
                }

                @Override
                public void mouseWheelMoved(MouseWheelEvent e) {
                    if (geometry == null) return;

                    // 计算缩放因子：滚轮向上滑动放大(1.1)，向下滑动缩小(0.9)
                    double zoomFactor = e.getWheelRotation() < 0 ? 1.1 : 0.9;

                    // 获取鼠标在屏幕上的位置
                    double mouseScreenX = e.getX();
                    double mouseScreenY = e.getY();

                    // 将鼠标屏幕坐标转换为几何坐标
                    double[] geomCoord = screenToGeom(mouseScreenX, mouseScreenY);
                    double mouseGeomX = geomCoord[0];
                    double mouseGeomY = geomCoord[1];

                    // 应用缩放
                    scale *= zoomFactor;

                    // 调整视图中心，使鼠标指向的几何点保持在屏幕相同位置
                    // 新视图中心 = 鼠标指向的几何点 - (鼠标屏幕偏移 / 新缩放比例)
                    double screenOffsetX = mouseScreenX - getWidth() / 2.0;
                    double screenOffsetY = mouseScreenY - getHeight() / 2.0;
                    viewCenterX = mouseGeomX - screenOffsetX / (autoScale * scale);
                    viewCenterY = mouseGeomY + screenOffsetY / (autoScale * scale); // Y轴翻转

                    repaint();
                }
            };

            addMouseListener(mouseAdapter);
            addMouseMotionListener(mouseAdapter);
            addMouseWheelListener(mouseAdapter);
        }

        /**
         * 将屏幕坐标转换为几何坐标
         */
        private double[] screenToGeom(double screenX, double screenY) {
            double centerX = getWidth() / 2.0;
            double centerY = getHeight() / 2.0;
            double geomX = viewCenterX + (screenX - centerX) / (autoScale * scale);
            double geomY = viewCenterY - (screenY - centerY) / (autoScale * scale); // Y轴翻转
            return new double[]{geomX, geomY};
        }

        public void setGeometry(Geometry geometry) {
            this.geometry = geometry;
            resetView();
            repaint();
        }

        private void resetView() {
            if (geometry != null) {
                Envelope envelope = geometry.getEnvelopeInternal();
                geomCenterX = (envelope.getMinX() + envelope.getMaxX()) / 2.0;
                geomCenterY = (envelope.getMinY() + envelope.getMaxY()) / 2.0;
                viewCenterX = geomCenterX;
                viewCenterY = geomCenterY;
            }
            scale = 1.0;
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            if (geometry == null) {
                g2d.setColor(Color.GRAY);
                g2d.setFont(new Font(Font.SANS_SERIF, Font.ITALIC, 14));
                String msg = "请输入WKT字符串并点击显示轮廓";
                FontMetrics fm = g2d.getFontMetrics();
                int x = (getWidth() - fm.stringWidth(msg)) / 2;
                int y = getHeight() / 2;
                g2d.drawString(msg, x, y);
                return;
            }

            // 计算几何图形的边界
            Envelope envelope = geometry.getEnvelopeInternal();
            double geomWidth = envelope.getWidth();
            double geomHeight = envelope.getHeight();

            if (geomWidth == 0 || geomHeight == 0) {
                return;
            }

            // 计算自动缩放比例以适应画布
            int padding = 50;
            double availableWidth = getWidth() - 2 * padding;
            double availableHeight = getHeight() - 2 * padding;

            double scaleX = availableWidth / geomWidth;
            double scaleY = availableHeight / geomHeight;
            autoScale = Math.min(scaleX, scaleY);

            // 保存变换
            AffineTransform oldTransform = g2d.getTransform();

            // 计算变换参数
            double centerX = getWidth() / 2.0;
            double centerY = getHeight() / 2.0;

            // 应用变换：
            // 1. 移动到画布中心
            // 2. 应用用户缩放
            // 3. Y轴翻转（因为屏幕Y轴向下，几何Y轴向上）
            // 4. 应用视图中心偏移
            g2d.translate(centerX, centerY);
            g2d.scale(autoScale * scale, -autoScale * scale);
            g2d.translate(-viewCenterX, -viewCenterY);

            // 绘制几何图形
            drawGeometry(g2d, geometry);

            // 恢复变换
            g2d.setTransform(oldTransform);
        }

        private void drawGeometry(Graphics2D g2d, Geometry geometry) {
            if (geometry instanceof Polygon) {
                drawPolygon(g2d, (Polygon) geometry);
            } else if (geometry instanceof MultiPolygon) {
                MultiPolygon multiPolygon = (MultiPolygon) geometry;
                for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                    Geometry geom = multiPolygon.getGeometryN(i);
                    if (geom instanceof Polygon) {
                        drawPolygon(g2d, (Polygon) geom);
                    }
                }
            }
        }

        private void drawPolygon(Graphics2D g2d, Polygon polygon) {
            // 绘制外轮廓填充
            Path2D exteriorPath = createPath(polygon.getExteriorRing().getCoordinates());
            g2d.setColor(new Color(255, 215, 0, 128)); // 半透明金色填充
            g2d.fill(exteriorPath);

            // 绘制内孔洞（用白色填充挖空）
            for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
                Path2D holePath = createPath(polygon.getInteriorRingN(i).getCoordinates());
                g2d.setColor(Color.WHITE);
                g2d.fill(holePath);
            }

            // 绘制外轮廓边框
            g2d.setColor(new Color(255, 140, 0)); // 深橙色边框
            g2d.setStroke(new BasicStroke(2.0f / (float) scale));
            g2d.draw(exteriorPath);

            // 绘制内孔洞边框
            g2d.setColor(new Color(200, 100, 0));
            g2d.setStroke(new BasicStroke(1.0f / (float) scale));
            for (int i = 0; i < polygon.getNumInteriorRing(); i++) {
                Path2D holePath = createPath(polygon.getInteriorRingN(i).getCoordinates());
                g2d.draw(holePath);
            }
        }

        private Path2D createPath(Coordinate[] coordinates) {
            Path2D path = new Path2D.Double();
            if (coordinates.length > 0) {
                path.moveTo(coordinates[0].x, coordinates[0].y);
                for (int i = 1; i < coordinates.length; i++) {
                    path.lineTo(coordinates[i].x, coordinates[i].y);
                }
                path.closePath();
            }
            return path;
        }
    }

    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            e.printStackTrace();
        }

        SwingUtilities.invokeLater(() -> {
            ShowGaussPointsGui gui = new ShowGaussPointsGui();
            gui.setVisible(true);
        });
    }
}
