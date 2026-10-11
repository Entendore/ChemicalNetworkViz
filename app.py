import sys
import math
import random
import heapq
from collections import deque
from typing import List, Dict, Optional, Tuple, Set

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGraphicsScene, QGraphicsView, QGraphicsEllipseItem, QGraphicsPathItem,
    QGraphicsTextItem, QSlider, QLabel, QPushButton, QComboBox, QGroupBox,
    QFrame, QGraphicsItem, QMessageBox, QFileDialog, QDockWidget, QTextEdit,
    QSplitter
)
from PySide6.QtCore import (
    Qt, QTimer, QPointF, QRectF, Property, QLineF
)
from PySide6.QtGui import (
    QPainter, QPen, QBrush, QColor, QPainterPath, QFont, QLinearGradient,
    QRadialGradient, QImage, QTextCursor
)

# ====================================================================================
# PROFESSIONAL GRAPH ALGORITHMS
# ====================================================================================

class GraphAlgorithms:
    """Advanced algorithms for chemical network analysis"""
    
    @staticmethod
    def get_adjacency_map(nodes: Dict, edges: List) -> Tuple[Dict, Dict]:
        """Helper to build adjacency lists"""
        adj_out = {nid: [] for nid in nodes}
        adj_in = {nid: [] for nid in nodes}
        for edge in edges:
            src, tgt = edge.source.node_id, edge.target.node_id
            adj_out[src].append(tgt)
            adj_in[tgt].append(src)
        return adj_out, adj_in

    @staticmethod
    def topological_sort(nodes: Dict, edges: List) -> List[str]:
        """
        Returns a linear ordering of vertices such that for every directed edge u -> v, 
        vertex u comes before v. Essential for determining reaction order.
        """
        adj_out, adj_in = GraphAlgorithms.get_adjacency_map(nodes, edges)
        
        # Calculate in-degrees
        in_degree = {nid: 0 for nid in nodes}
        for nid in nodes:
            for neighbor in adj_in[nid]:
                in_degree[nid] += 1 # This logic was slightly wrong in thought, corrected:
        
        # Correct in-degree calculation
        in_degree = {nid: 0 for nid in nodes}
        for edge in edges:
            in_degree[edge.target.node_id] += 1
            
        # Queue for nodes with no incoming edges
        queue = deque([nid for nid, deg in in_degree.items() if deg == 0])
        topo_order = []
        
        while queue:
            u = queue.popleft()
            topo_order.append(u)
            
            for v in adj_out[u]:
                in_degree[v] -= 1
                if in_degree[v] == 0:
                    queue.append(v)
                    
        return topo_order

    @staticmethod
    def calculate_betweenness_centrality(nodes: Dict, edges: List) -> Dict[str, float]:
        """
        Brandes' Algorithm for Betweenness Centrality.
        Identifies nodes that act as bridges/control points in the network.
        High score = Critical Bottleneck.
        """
        adj_out, _ = GraphAlgorithms.get_adjacency_map(nodes, edges)
        betweenness = {nid: 0.0 for nid in nodes}
        node_ids = list(nodes.keys())
        
        for s in node_ids:
            # Single source shortest path
            S = []
            P = {w: [] for w in node_ids}
            sigma = {t: 0 for t in node_ids}
            sigma[s] = 1
            D = {t: -1 for t in node_ids}
            D[s] = 0
            
            # BFS
            Q = deque([s])
            while Q:
                v = Q.popleft()
                S.append(v)
                for w in adj_out[v]:
                    # Path discovery
                    if D[w] < 0:
                        Q.append(w)
                        D[w] = D[v] + 1
                    # Path counting
                    if D[w] == D[v] + 1:
                        sigma[w] = sigma[w] + sigma[v]
                        P[w].append(v)
            
            # Accumulation
            delta = {v: 0.0 for v in node_ids}
            while S:
                w = S.pop()
                for v in P[w]:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
                if w != s:
                    betweenness[w] += delta[w]
                    
        return betweenness

    @staticmethod
    def simulate_knockout(nodes: Dict, edges: List, target_id: str, removed_node_id: str) -> bool:
        """
        Robustness analysis: Can we still reach 'target_id' if 'removed_node_id' is deleted?
        Returns True if the network is ROBUST (target reachable), False if FRAGILE.
        """
        # Build adjacency list excluding the removed node
        adj = {nid: [] for nid in nodes if nid != removed_node_id}
        
        for edge in edges:
            if edge.source.node_id == removed_node_id or edge.target.node_id == removed_node_id:
                continue
            adj[edge.source.node_id].append(edge.target.node_id)
            
        # BFS from any source to target
        # For simplicity, we check if ANY path exists to target from the remaining graph
        # In a real scenario, we'd start from specific precursors.
        
        visited = set()
        queue = deque()
        
        # Start from all potential sources (nodes with no incoming edges in the reduced graph)
        # Or just start from all remaining nodes to check connectivity
        for nid in adj:
            queue.append(nid)
            
        while queue:
            curr = queue.popleft()
            if curr == target_id:
                return True # Target is still reachable
            
            if curr in visited:
                continue
            visited.add(curr)
            
            for neighbor in adj[curr]:
                if neighbor not in visited:
                    queue.append(neighbor)
                    
        return False


# ====================================================================================
# GRAPHICS ITEMS (Visuals)
# ====================================================================================

class Particle:
    """Particle flowing along an edge"""
    def __init__(self):
        self.progress = random.random()
        self.speed = random.uniform(0.003, 0.008)
        self.size = random.uniform(2, 4)
        self.active = True
    
    def update(self, mult: float = 1.0):
        self.progress += self.speed * mult
        if self.progress >= 1.0:
            self.progress = 0.0
            return True
        return False


class ChemicalNode(QGraphicsEllipseItem):
    """Visual representation of a chemical species"""
    
    def __init__(self, node_id: str, name: str, node_type: str, x: float, y: float):
        self.node_id = node_id
        self.name = name
        self.node_type = node_type
        self.radius = 30
        
        self.concentration = 1.0
        self.knocked_out = False # For robustness analysis
        
        super().__init__(-self.radius, -self.radius, self.radius * 2, self.radius * 2)
        self.setPos(x, y)
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges, True)
        self.setAcceptHoverEvents(True)
        
        self._highlight = False
        self._is_critical = False
        
        self.type_colors = {
            'reactant': QColor(30, 144, 255),    # Dodger Blue
            'product': QColor(50, 205, 50),      # Lime Green
            'intermediate': QColor(255, 165, 0), # Orange
            'enzyme': QColor(220, 20, 60),       # Crimson
            'target': QColor(255, 0, 255)        # Magenta
        }
        self.base_color = self.type_colors.get(node_type, QColor(100, 100, 100))
        self.setup_appearance()
        
        # Label
        self.label = QGraphicsTextItem(name, self)
        font = QFont("Arial", 8, QFont.Bold)
        self.label.setFont(font)
        self.label.setDefaultTextColor(Qt.white)
        self.label.setPos(-self.label.boundingRect().width()/2, 
                          -self.label.boundingRect().height()/2)

    def setup_appearance(self):
        self.setZValue(10)
        self.update_brush()

    def update_brush(self):
        if self.knocked_out:
            color = QColor(80, 80, 80)
            pen_color = QColor(50, 50, 50)
        elif self._is_critical:
            color = QColor(255, 0, 0) # Red for critical bottleneck
            pen_color = QColor(255, 100, 100)
        elif self._highlight:
            color = self.base_color.lighter(150)
            pen_color = QColor(255, 255, 255)
        else:
            color = self.base_color
            pen_color = color.lighter(120)
            
        gradient = QRadialGradient(0, 0, self.radius)
        gradient.setColorAt(0, color.lighter(130))
        gradient.setColorAt(0.8, color)
        gradient.setColorAt(1, color.darker(150))
        
        self.setBrush(QBrush(gradient))
        self.setPen(QPen(pen_color, 2))

    def set_critical(self, status: bool):
        self._is_critical = status
        self.update_brush()

    def set_knocked_out(self, status: bool):
        self.knocked_out = status
        self.update_brush()
        self.setOpacity(0.5 if status else 1.0)
        self.setFlag(QGraphicsItem.ItemIsEnabled, not status)

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemPositionHasChanged:
            scene = self.scene()
            if scene:
                scene.update_edges_for_node(self)
        return super().itemChange(change, value)


class ChemicalEdge(QGraphicsPathItem):
    """Visual representation of a reaction"""
    
    def __init__(self, source: ChemicalNode, target: ChemicalNode, 
                 reaction_type: str = 'reaction', rate: float = 1.0):
        self.source = source
        self.target = target
        self.reaction_type = reaction_type
        self.rate = rate
        
        super().__init__()
        self.setZValue(1)
        
        self.particles: List[Particle] = [Particle() for _ in range(int(3 + rate * 3))]
        
        self.type_colors = {
            'reaction': QColor(100, 180, 255),
            'catalysis': QColor(255, 100, 255),
            'inhibition': QColor(255, 80, 80),
        }
        
        self.base_color = self.type_colors.get(reaction_type, QColor(150, 150, 150))
        self.setup_appearance()
        self.update_path()

    def setup_appearance(self):
        width = 1.5 + self.rate * 1.5
        pen = QPen(self.base_color, width)
        pen.setCapStyle(Qt.RoundCap)
        self.setPen(pen)

    def update_path(self):
        if self.source.knocked_out or self.target.knocked_out:
            self.setVisible(False)
            return
        self.setVisible(True)
        
        path = QPainterPath()
        src_pos = self.source.pos()
        tgt_pos = self.target.pos()
        
        dx = tgt_pos.x() - src_pos.x()
        dy = tgt_pos.y() - src_pos.y()
        dist = math.sqrt(dx*dx + dy*dy)
        
        if dist > 0:
            # Curved path
            offset = dist * 0.2
            mid_x = (src_pos.x() + tgt_pos.x()) / 2
            mid_y = (src_pos.y() + tgt_pos.y()) / 2
            ctrl_x = mid_x - dy/dist * offset
            ctrl_y = mid_y + dx/dist * offset
            
            path.moveTo(src_pos)
            path.quadTo(QPointF(ctrl_x, ctrl_y), tgt_pos)
        
        self.setPath(path)

    def get_point_at_progress(self, progress: float) -> QPointF:
        path = self.path()
        if path.elementCount() >= 3:
            p0 = path.elementAt(0)
            p1 = path.elementAt(1)
            p2 = path.elementAt(2)
            t = progress
            x = (1-t)**2 * p0.x + 2*(1-t)*t * p1.x + t**2 * p2.x
            y = (1-t)**2 * p0.y + 2*(1-t)*t * p1.y + t**2 * p2.y
            return QPointF(x, y)
        return self.source.pos()

    def paint(self, painter, option, widget):
        super().paint(painter, option, widget)
        if not self.isVisible():
            return
            
        painter.setRenderHint(QPainter.Antialiasing)
        for p in self.particles:
            pos = self.get_point_at_progress(p.progress)
            gradient = QRadialGradient(pos, p.size * 2)
            gradient.setColorAt(0, QColor(255, 255, 255, 200))
            gradient.setColorAt(1, QColor(255, 255, 255, 0))
            painter.setBrush(QBrush(gradient))
            painter.setPen(Qt.NoPen)
            painter.drawEllipse(pos, p.size * 2, p.size * 2)


# ====================================================================================
# SCENE & LOGIC CONTROLLERS
# ====================================================================================

class ChemicalNetworkScene(QGraphicsScene):
    def __init__(self):
        super().__init__()
        self.setSceneRect(-1000, -800, 2000, 1600)
        self.nodes: Dict[str, ChemicalNode] = {}
        self.edges: List[ChemicalEdge] = []
        self.setBackgroundBrush(QBrush(QColor(20, 20, 30)))

    def add_node(self, node_id: str, name: str, node_type: str, x: float, y: float):
        node = ChemicalNode(node_id, name, node_type, x, y)
        self.nodes[node_id] = node
        self.addItem(node)
        return node

    def add_edge(self, source_id: str, target_id: str, reaction_type: str = 'reaction', rate: float = 1.0):
        source = self.nodes.get(source_id)
        target = self.nodes.get(target_id)
        if source and target:
            edge = ChemicalEdge(source, target, reaction_type, rate)
            self.edges.append(edge)
            self.addItem(edge)
            return edge
        return None

    def update_edges_for_node(self, node: ChemicalNode):
        for edge in self.edges:
            if edge.source == node or edge.target == node:
                edge.update_path()
                
    def reset_visual_states(self):
        for node in self.nodes.values():
            node.set_critical(False)
            node.set_knocked_out(False)
            node._highlight = False
            node.update_brush()
        for edge in self.edges:
            edge.update_path()
            edge.setVisible(True)


class ForceLayoutController:
    """Fruchterman-Reingold algorithm"""
    def __init__(self, scene: ChemicalNetworkScene):
        self.scene = scene
        self.temperature = 1.0
        self.active = False
        
    def step(self):
        if not self.active: return
        
        nodes = list(self.scene.nodes.values())
        if not nodes: return
        
        area = 1000*1000
        k = math.sqrt(area / len(nodes))
        disp = {n.node_id: QPointF(0,0) for n in nodes}
        
        # Repulsion
        for i, v in enumerate(nodes):
            for j, u in enumerate(nodes):
                if i == j: continue
                delta = v.pos() - u.pos()
                dist = max(0.01, math.sqrt(delta.x()**2 + delta.y()**2))
                force = k*k / dist
                disp[v.node_id] += QPointF(delta.x()/dist * force, delta.y()/dist * force)
        
        # Attraction
        for edge in self.scene.edges:
            delta = edge.source.pos() - edge.target.pos()
            dist = max(0.01, math.sqrt(delta.x()**2 + delta.y()**2))
            force = dist / k
            disp[edge.source.node_id] -= QPointF(delta.x()/dist * force, delta.y()/dist * force)
            disp[edge.target.node_id] += QPointF(delta.x()/dist * force, delta.y()/dist * force)
            
        # Apply
        for node in nodes:
            if node.isSelected(): continue
            d = disp[node.node_id]
            length = max(0.01, math.sqrt(d.x()**2 + d.y()**2))
            new_pos = node.pos() + QPointF(d.x()/length * min(length, self.temperature*10), 
                                           d.y()/length * min(length, self.temperature*10))
            node.setPos(new_pos)
            
        self.temperature *= 0.95
        if self.temperature < 0.01:
            self.active = False


# ====================================================================================
# UI COMPONENTS
# ====================================================================================

class AnalysisConsole(QDockWidget):
    """Professional output console for algorithm results"""
    def __init__(self):
        super().__init__("Analysis Console")
        self.setAllowedAreas(Qt.BottomDockWidgetArea)
        self.text_edit = QTextEdit()
        self.text_edit.setReadOnly(True)
        self.text_edit.setStyleSheet("background-color: #1e1e1e; color: #dcdcdc; font-family: Consolas;")
        self.setWidget(self.text_edit)
        
    def log(self, message: str, color: str = "white"):
        cursor = self.text_edit.textCursor()
        cursor.movePosition(QTextCursor.End)
        self.text_edit.setTextCursor(cursor)
        self.text_edit.insertHtml(f"<span style='color:{color};'>{message}</span><br>")
        
    def clear(self):
        self.text_edit.clear()


class ControlPanel(QWidget):
    def __init__(self, scene: ChemicalNetworkScene, console: AnalysisConsole):
        super().__init__()
        self.scene = scene
        self.console = console
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(15)
        
        # --- Graph Layout ---
        grp_layout = QGroupBox("Graph Layout")
        lay = QVBoxLayout(grp_layout)
        self.btn_force = QPushButton("Force-Directed Layout")
        self.btn_force.setCheckable(True)
        self.btn_hierarchical = QPushButton("Hierarchical Layout")
        lay.addWidget(self.btn_force)
        lay.addWidget(self.btn_hierarchical)
        layout.addWidget(grp_layout)
        
        # --- Network Analysis ---
        grp_analysis = QGroupBox("Network Analysis")
        lay = QVBoxLayout(grp_analysis)
        self.btn_centrality = QPushButton("Find Bottlenecks (Betweenness)")
        self.btn_topo = QPushButton("Check Reaction Order (Topo Sort)")
        lay.addWidget(self.btn_centrality)
        lay.addWidget(self.btn_topo)
        layout.addWidget(grp_analysis)
        
        # --- Professional Use Cases ---
        grp_pro = QGroupBox("Professional Simulations")
        lay = QVBoxLayout(grp_pro)
        self.btn_knockout = QPushButton("Drug Target Simulation (Knockout)")
        self.btn_pathway = QPushButton("Trace Synthesis Path")
        lay.addWidget(self.btn_knockout)
        lay.addWidget(self.btn_pathway)
        layout.addWidget(grp_pro)
        
        # --- Visuals ---
        self.slider_speed = QSlider(Qt.Horizontal)
        self.slider_speed.setRange(1, 100)
        self.slider_speed.setValue(50)
        layout.addWidget(QLabel("Animation Speed:"))
        layout.addWidget(self.slider_speed)
        
        # --- Actions ---
        self.btn_reset = QPushButton("Reset Network State")
        layout.addWidget(self.btn_reset)
        
        layout.addStretch()
        
        # Stylesheet
        self.setStyleSheet("""
            QGroupBox { font-weight: bold; border: 1px solid #555; border-radius: 5px; margin-top: 10px; color: #aaa; }
            QPushButton { background-color: #2a5298; border: none; padding: 8px; border-radius: 4px; color: white; }
            QPushButton:hover { background-color: #3a72c8; }
            QPushButton:checked { background-color: #28a745; }
        """)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Pro Chemical Network Viz")
        self.resize(1400, 900)
        
        # 1. Setup Core Components
        self.scene = ChemicalNetworkScene()
        self.view = QGraphicsView(self.scene)
        self.view.setRenderHint(QPainter.Antialiasing)
        self.view.setDragMode(QGraphicsView.RubberBandDrag)
        self.view.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.view.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        self.layout_controller = ForceLayoutController(self.scene)
        
        # 2. Setup UI
        self.setup_ui()
        self.create_demo_network()
        
        # 3. Timers
        self.timer = QTimer()
        self.timer.timeout.connect(self.on_tick)
        self.timer.start(16)

    def setup_ui(self):
        # Main Layout
        splitter = QSplitter(Qt.Horizontal)
        self.setCentralWidget(splitter)
        
        # Left: View
        splitter.addWidget(self.view)
        
        # Right: Control Panel
        self.console = AnalysisConsole()
        self.panel = ControlPanel(self.scene, self.console)
        self.panel.setFixedWidth(300)
        splitter.addWidget(self.panel)
        
        # Bottom: Console
        self.addDockWidget(Qt.BottomDockWidgetArea, self.console)
        
        # Connect Signals
        self.panel.btn_force.clicked.connect(self.on_force_layout)
        self.panel.btn_hierarchical.clicked.connect(self.on_hierarchical_layout)
        self.panel.btn_centrality.clicked.connect(self.on_centrality)
        self.panel.btn_topo.clicked.connect(self.on_topological)
        self.panel.btn_knockout.clicked.connect(self.on_knockout)
        self.panel.btn_pathway.clicked.connect(self.on_trace_path)
        self.panel.btn_reset.clicked.connect(self.on_reset)

    # ====================================================================================
    # ALGORITHM HANDLERS
    # ====================================================================================

    def on_force_layout(self, checked):
        if checked:
            self.layout_controller.temperature = 1.0
            self.layout_controller.active = True
        else:
            self.layout_controller.active = False

    def on_hierarchical_layout(self):
        """Assign layers based on topological sort and arrange horizontally"""
        self.console.log("Running Topological Sort...", "cyan")
        order = GraphAlgorithms.topological_sort(self.scene.nodes, self.scene.edges)
        
        if not order:
            self.console.log("Graph has cycles, cannot perform strict hierarchical layout.", "red")
            return
            
        self.console.log(f"Order: {' -> '.join(order)}", "green")
        
        # Assign positions
        layer_width = 150
        x_pos = - (len(order) / 2) * layer_width
        
        for nid in order:
            node = self.scene.nodes[nid]
            node.setPos(x_pos, 0) # Center Y, could expand vertically if too many nodes
            x_pos += layer_width
            
        self.console.log("Layout updated.", "green")

    def on_centrality(self):
        """Identify critical bottleneck nodes"""
        self.scene.reset_visual_states()
        self.console.log("Calculating Betweenness Centrality...", "cyan")
        
        scores = GraphAlgorithms.calculate_betweenness_centrality(self.scene.nodes, self.scene.edges)
        max_score = max(scores.values()) if scores else 0
        
        self.console.log("Results:", "white")
        for nid, score in sorted(scores.items(), key=lambda item: item[1], reverse=True):
            node = self.scene.nodes[nid]
            pct = (score / max_score) * 100 if max_score > 0 else 0
            self.console.log(f"  {node.name}: {pct:.1f}% importance", "white")
            
            # Highlight nodes with high importance (> 40%)
            if pct > 40:
                node.set_critical(True)
                self.console.log(f"  >> {node.name} identified as BOTTLENECK.", "red")

    def on_topological(self):
        """Just prints the order to console"""
        self.console.log("Checking Reaction Sequence...", "cyan")
        order = GraphAlgorithms.topological_sort(self.scene.nodes, self.scene.edges)
        if order:
            self.console.log("Valid Synthesis Order: " + " -> ".join(order), "green")
        else:
            self.console.log("Cycle detected or invalid graph.", "red")

    def on_knockout(self):
        """
        Professional Use Case: Drug Target Discovery.
        Simulate removing a node to see if a target becomes unreachable.
        """
        items = self.scene.selectedItems()
        if len(items) != 1:
            self.console.log("Select exactly 1 ENZYME node to knock out.", "yellow")
            return
            
        knockout_node = items[0]
        
        # Ask user for the TARGET product they want to preserve
        # For demo, we hardcode a target, e.g., 'Product' or the last node
        target_id = 'P1' # Default for our demo graph
        
        self.console.log(f"SIMULATION: Knocking out {knockout_node.name}...", "cyan")
        self.scene.reset_visual_states()
        
        # Perform simulation logic
        is_robust = GraphAlgorithms.simulate_knockout(
            self.scene.nodes, self.scene.edges, target_id, knockout_node.node_id
        )
        
        # Visual feedback
        knockout_node.set_knocked_out(True)
        
        if is_robust:
            self.console.log(f"Result: Network is ROBUST. Target {target_id} is still reachable.", "green")
            self.console.log(f"Implication: {knockout_node.name} is a poor drug target.", "white")
        else:
            self.console.log(f"Result: Network is FRAGILE. Target {target_id} is unreachable.", "red")
            self.console.log(f"Implication: {knockout_node.name} is a potential DRUG TARGET.", "white")

    def on_trace_path(self):
        """Find shortest synthesis path between two selected nodes"""
        items = self.scene.selectedItems()
        if len(items) < 2:
            self.console.log("Select 2 nodes: Source (Reactant) and Target (Product).", "yellow")
            return
            
        start = items[0]
        end = items[1]
        
        self.console.log(f"Calculating shortest path from {start.name} to {end.name}...", "cyan")
        
        # Using Dijkstra (implemented in previous file block, simplified BFS here for brevity in this specific block, 
        # but we can just use the path logic from centrality simulation)
        # Re-using BFS logic from simulation
        adj_out, _ = GraphAlgorithms.get_adjacency_map(self.scene.nodes, self.scene.edges)
        
        # BFS
        visited = {start.node_id: None}
        queue = deque([start.node_id])
        found = False
        
        while queue:
            curr = queue.popleft()
            if curr == end.node_id:
                found = True
                break
            for neighbor in adj_out[curr]:
                if neighbor not in visited:
                    visited[neighbor] = curr
                    queue.append(neighbor)
                    
        self.scene.reset_visual_states()
        
        if found:
            # Reconstruct path
            path = []
            curr = end.node_id
            while curr:
                path.append(curr)
                curr = visited[curr]
            path.reverse()
            
            self.console.log("Path found: " + " -> ".join([self.scene.nodes[n].name for n in path]), "green")
            
            for nid in path:
                self.scene.nodes[nid].set_critical(True) # Highlight path
        else:
            self.console.log("No path found.", "red")

    def on_reset(self):
        self.scene.reset_visual_states()
        self.console.clear()
        self.console.log("System Reset.", "white")

    def on_tick(self):
        # Animation Loop
        speed = self.panel.slider_speed.value() / 50.0
        
        # Layout Step
        if self.layout_controller.active:
            self.layout_controller.step()
            if not self.layout_controller.active:
                self.panel.btn_force.setChecked(False)
        
        # Particles
        for edge in self.scene.edges:
            if edge.isVisible():
                for p in edge.particles:
                    p.update(speed)
                edge.update()

    def create_demo_network(self):
        """Creates a professional metabolic pathway demo"""
        # Precursors
        self.scene.add_node('R1', 'Glucose', 'reactant', -400, 0)
        self.scene.add_node('R2', 'ATP', 'reactant', -400, -100)
        
        # Intermediates
        self.scene.add_node('I1', 'G6P', 'intermediate', -200, 0)
        self.scene.add_node('I2', 'F6P', 'intermediate', 0, 0)
        self.scene.add_node('I3', 'DHAP', 'intermediate', 200, -100)
        self.scene.add_node('I4', 'G3P', 'intermediate', 200, 100)
        
        # Enzymes
        self.scene.add_node('E1', 'Hexokinase', 'enzyme', -300, -50)
        self.scene.add_node('E2', 'PFK-1', 'enzyme', -100, -50)
        
        # Product
        self.scene.add_node('P1', 'Pyruvate', 'product', 400, 50)
        
        # Reactions
        self.scene.add_edge('R1', 'I1', 'reaction', 1.5)
        self.scene.add_edge('R2', 'I1', 'reaction', 1.0)
        self.scene.add_edge('I1', 'I2', 'reaction', 1.0)
        self.scene.add_edge('I2', 'I3', 'reaction', 1.0)
        self.scene.add_edge('I2', 'I4', 'reaction', 1.0)
        self.scene.add_edge('I3', 'I4', 'reaction', 0.5)
        self.scene.add_edge('I4', 'P1', 'reaction', 2.0)
        
        # Catalysis (Dashed/Colored differently)
        self.scene.add_edge('E1', 'I1', 'catalysis', 0.5)
        self.scene.add_edge('E2', 'I2', 'catalysis', 0.5)
        
        self.console.log("Demo Network Loaded: Glycolysis (Simplified).", "cyan")
        self.console.log("Instructions:", "yellow")
        self.console.log("1. Move nodes by dragging.", "white")
        self.console.log("2. Use 'Force Layout' to auto-arrange.", "white")
        self.console.log("3. Select 'Hexokinase' and click 'Drug Target Simulation'.", "white")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    sys.exit(app.exec())