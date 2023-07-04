import tester.*;
import javalib.impworld.*;
import java.awt.Color;
import javalib.worldimages.*;
import java.util.*;
import java.util.function.Predicate;

/*
 * MAZE INSTRUCTIONS:
 * 
 * B - runs breadth first search of maze
 * D - runs depth first search of maze
 * M - allows user to use arrow keys to solve maze manually (EXTRA CREDIT)
 * up arrow - moves current vertex up
 * down arrow - moves current vertex down
 * left arrow - moves current vertex left
 * right arrow - moves current vertex right
 * R - resets maze without exiting program (EXTRA CREDIT) 
 * */


// function object to compare the weight of edges 
// this will be used to implement Kruskal's algorithm
// for maze construction

class SortEdges implements Comparator<Edge> {

  public int compare(Edge e1, Edge e2) {
    if (e1.weight < e2.weight) {
      return 1;
    } else if (e1.weight > e2.weight) {
      return -1;
    } else {
      return 0;
    }
  }
}

//function object for the Predicate for Strings

class StringCheck implements Predicate<String> {
  String s1;

  // the constructor

  StringCheck(String s1) {
    this.s1 = s1;
  }

  // checks whether this and a given string are equal

  public boolean test(String t) {
    return s1.equals(t);
  }
}


// represents a vertex in a graph

class Vertex {
  int x;
  int y;
  boolean visited;
  boolean solution;
  ArrayList<Edge> neighborEdges;

  Vertex(int x, int y) {
    this.x = x;
    this.y = y;
    this.neighborEdges = new ArrayList<Edge>();
    this.visited = false;
    this.solution = false;
  }

  // used for identifying each vertex
  
  int identifyVertex() {
    return 1000 * y + x;
  }
}

// represents an edge in a graph

class Edge {
  Vertex from;
  Vertex to;
  int weight;

  Edge(Vertex from, Vertex to, int weight) {
    this.from = from;
    this.to = to;
    this.weight = weight;
  }
}


// represents a deque 

class Deque<T> {
  Sentinel<T> header;

  Deque() {
    this.header = new Sentinel<T>();
  }

  Deque(Sentinel<T> header) {
    this.header = header;
  }

  // gets the size of a deque
  int size() {
    return this.header.size();
  }

  // adds a node to the head of a deque
  T addAtHead(T t) {
    return header.addAtHead(t);
  }

  // adds a node to the tail of a deque
  T addAtTail(T t) {
    return header.addAtTail(t);
  }
  
  // removes a node from a deque
  void removeNode(ANode<T> aNode) {
    if (!aNode.equals(header)) {
      header.removeNodeFirst(aNode);
    }
  }

  // removes a node from the head
  T removeFromHead() {
    return header.removeFromHead();
  }

  // removes a node from the tail
  T removeFromTail() {
    return header.removeFromTail();
  }

  // gets a node that satisfies the predicate
  
  ANode<T> find(Predicate<T> p) {
    return header.find(p);
  }

}


// represents either a Node or Sientel

abstract class ANode<T> {
  ANode<T> next;
  ANode<T> prev;

  // returns the total number of nodes in a deque
  
  int size(ANode<T> s) {
    if (this.next.equals(s)) {
      return 1;
    }
    else {
      return 1 + this.next.size(s);
    }
  }

  // removes a node from a deque
  T remove() {
    this.prev.next = this.next;
    this.next.prev = this.prev;
    return this.getData();
  }

  // finds a node that follows the predicate
  
  ANode<T> find(Predicate<T> p) {
    return this;
  }

  // helper to remove a node
  
  void removeNode(ANode<T> aNode) {
    if (this.equals(aNode)) {
      this.remove();
    }
    else {
      this.next.removeNode(aNode);
    }
  }

  // returns the node data
  
  T getData() {
    return null;
  }

  // checks whether this node has a next value
  
  boolean hasNextVal() {
    return !(next.isSentinel());
  }

  // checks if an ANode is a Sentinel
  boolean isSentinel() {
    return false;
  }
}

// represents a Sentinel

class Sentinel<T> extends ANode<T> {

  Sentinel() {
    this.next = this;
    this.prev = this;
  }

  // gets the size of a sentinel
  
  int size() {
    if (this.next.equals(this)) {
      return 0;
    }
    return this.next.size(this);
  }

  // adds a node with data t to the head
  
  T addAtHead(T t) {
    new Node<T>(t, this.next, this);
    return this.next.getData();
  }

  // adds a node with data t to the tail
  
  T addAtTail(T t) {
    new Node<T>(t, this, this.prev);
    return this.next.getData();
  }

  // removes a node from the head of a deque
  
  T removeFromHead() {
    return this.next.remove();
  }

  // removes a node from the tail of a deque
  
  T removeFromTail() {
    return this.prev.remove();
  }

  // finds the next Anode that follows a given predicate
  
  ANode<T> find(Predicate<T> p) {
    return this.next.find(p);
  }


  // removes the first node
  
  void removeNodeFirst(ANode<T> aNode) {
    this.next.removeNode(aNode);
  }


  // helps removes a node
  
  void removeNode(ANode<T> aNode) {
    return;
  }

  // checks if this node is a sentinel
  
  boolean isSentinel() {
    return true;
  }
}

// represents a Node

class Node<T> extends ANode<T> {
  T data;

  Node(T data) {
    this.data = data;
    this.next = null;
    this.prev = null;
  }

  Node(T data, ANode<T> next, ANode<T> prev) {
    if ((next == null) || (prev == null)) {
      throw new IllegalArgumentException("Cannot accept null node");
    }
    this.data = data;
    this.next = next;
    this.prev = prev;
    prev.next = this;
    next.prev = this;
  }

  // finds the first node that follows the condition

  ANode<T> find(Predicate<T> pred) {
    if (pred.test(this.data)) {
      return this;
    } else {
      return this.next.find(pred);
    }
  }


  // returns the data for this node
  
  T getData() {
    return this.data;
  }
}


// represents a queue 

class Queue<T> {
  // the list of our queue
  Deque<T> list;

  Queue() {
    this.list = new Deque<T>();
  }

  Queue(ArrayList<T> ts) {
    list = new Deque<T>();
    for (T t : ts) {
      list.addAtTail(t);
    }
  }

  // adds an item to the tail of the list
  void add(T item) {
    list.addAtTail(item);
  }

  // determines if our list are empty
  boolean isEmpty() {
    return list.size() == 0;
  }

  // removes and returns the tail of the list
  T remove() {
    return list.removeFromTail();
  }

  void empty() {
    list = new Deque<T>();
  }

}

// represents a search

abstract class Search {
  
  HashMap<Integer, Vertex> cameFromEdge;

  // gets the solution 
  
  void reconstruct(HashMap<Integer, Vertex> h, Vertex next) {
    while (h.containsKey(next.identifyVertex())) {
      next.solution = true;
      next = h.get(next.identifyVertex());
    }
  }
}

// represents a breadth first search and iterates through

class BFS extends Search {
  // our current worklist
  Queue<Vertex> worklist;

  BFS(ArrayList<Vertex> list) {
    this.worklist = new Queue<Vertex>();
    worklist.add(list.get(0));
    list.get(0).visited = true;
    cameFromEdge = new HashMap<Integer, Vertex>();
  }

  // determines if we have another task in our worklist
  public boolean hasNext() {
    return !worklist.isEmpty();
  }

  // goes through breadth first algorithm
  
  public Queue<Vertex> next() {
    Vertex u = worklist.remove();
    for (Edge e : u.neighborEdges) {
      if (!e.to.visited) {
        cameFromEdge.put(e.to.identifyVertex(), e.from);
        if (e.to.x == 50 - 1 && e.to.y == 30 - 1) {
          reconstruct(cameFromEdge, e.to);
          worklist = new Queue<Vertex>();
        }
        else {
          e.to.visited = true;
          worklist.add(e.to);
        }
      }
    }
    return worklist;
  }
}

// represents a depth first search and iterates through

class DFS extends Search {
  // our current worklist
  Stack<Vertex> worklist;

  DFS(ArrayList<Vertex> list) {
    this.worklist = new Stack<Vertex>();
    worklist.push(list.get(0));
    list.get(0).visited = true;
    cameFromEdge = new HashMap<Integer, Vertex>();
  }

  // checks if worklist is not empty
  
  public boolean hasNext() {
    return !worklist.isEmpty();
  }

  // goes through depth first algorithm
  
  public Stack<Vertex> next() {
    Vertex u = worklist.pop();
    for (Edge e : u.neighborEdges) {
      if (!e.to.visited) {
        cameFromEdge.put(e.to.identifyVertex(), e.from);
        if (e.to.x == 50 - 1 && e.to.y == 30 - 1) {
          reconstruct(cameFromEdge, e.to);
          worklist = new Stack<Vertex>();
        }
        else {
          worklist.push(u);
          e.to.visited = true;
          worklist.push(e.to);
          break;
        }
      }
    }
    return worklist;
  }
}


// EXTRA CREDIT - added feature for player to move in maze 
// with keyboard 

// class for iterating maze manually

class Manual extends Search {
  Vertex vertex;
  boolean done;

  Manual(ArrayList<Vertex> list) {
    vertex = list.get(0);
    cameFromEdge = new HashMap<Integer, Vertex>();
    done = false;
  }

  // checks if the maze is solved
  
  public boolean isSolved() {
    return !done;
  }

  // moves to next vertex in maze if not done
  
  public Vertex proceed(boolean b, Edge edge) {
    if (b) {
      vertex.visited = true;
      vertex.solution = false;
      if (!edge.to.visited) {
        cameFromEdge.put(edge.to.identifyVertex(), edge.from);
      }
      if (edge.to.x == 50 - 1 && edge.to.y == 30 - 1) {
        reconstruct(cameFromEdge, edge.to);
      }
      else {
        vertex = edge.to;
        vertex.solution = true;
      }
    }
    return vertex;
  }

  // moves the vertex left 
  
  public Vertex goLeft() {
    for (Edge e : vertex.neighborEdges) {
      this.proceed(e.to.x == vertex.x - 1, e);
    }
    return vertex;
  }

  // moves the vertex right 
  
  public Vertex goRight() {
    for (Edge e : vertex.neighborEdges) {
      this.proceed(e.to.x == vertex.x + 1, e);
    }
    return vertex;
  }

  // moves the vertex down  
  
  public Vertex goDown() {
    for (Edge e : vertex.neighborEdges) {
      this.proceed(e.to.y == vertex.y + 1, e);
    }
    return vertex;
  }

  public Vertex goUp() {
    for (Edge e : vertex.neighborEdges) {
      this.proceed(e.to.y == vertex.y - 1, e);
    }
    return vertex;
  }
}


// represents a Maze

class Maze extends World {
  int width = 50;
  int height = 30;
  int scale = 10;
  
  boolean bfs;
  boolean dfs;
  boolean manual;
  BFS b;
  DFS d;
  Manual m;
  // all of our vertices
  ArrayList<Vertex> vertices;
  // all of the walls
  ArrayList<Edge> walls;

  // constructs a maze
  Maze() {
    this.initMaze();
  }

  // sets up our maze by constructing the vertices, walls and initializing
  // our variables

  void initMaze() {
    ArrayList<ArrayList<Vertex>> v = makeGraph();
    ArrayList<Edge> edgeList = createEdges(v);
    v = kruskalAlgorithm(v);
    walls = createWalls(v, edgeList);
    vertices = new ArrayList<Vertex>();
    for (ArrayList<Vertex> vList : v) {
      for (Vertex vertex : vList) {
        vertices.add(vertex);
      }
    }
    bfs = false;
    dfs = false;
    manual = false;
  }

  // makes the walls for our maze by checking the upright walls

  ArrayList<Edge> createWalls(ArrayList<ArrayList<Vertex>> v, ArrayList<Edge> allEdge) {
    ArrayList<Edge> w = new ArrayList<Edge>();
    for (Edge e : allEdge) {
      boolean valid = true;
      for (ArrayList<Vertex> l : v) {
        for (Vertex vertex : l) {
          for (Edge e2 : vertex.neighborEdges) {
            if (e.equals(e2) || (e.to == e2.from && e.from == e2.to)) {
              valid = false;
            }
          }
        }
      }
      if (valid) {
        w.add(e);
      }
    }
    return w;
  }

  // returns an arraylist of every edge in a given graph
  // goes in the graph to get all edges and returns to user

  ArrayList<Edge> createEdges(ArrayList<ArrayList<Vertex>> graph) {
    ArrayList<Edge> allEdge = new ArrayList<Edge>();
    for (ArrayList<Vertex> verts : graph) {
      for (Vertex vertex : verts) {
        for (Edge edge : vertex.neighborEdges) {
          allEdge.add(edge);
        }
      }
    }
    return allEdge;
  }

  // creates aNode 2D grid of vertices which has edges of random weights

  ArrayList<ArrayList<Vertex>> makeGraph() {
    ArrayList<ArrayList<Vertex>> vertices = new ArrayList<ArrayList<Vertex>>();
    for (int x = 0; x < width; x++) {
      ArrayList<Vertex> temp = new ArrayList<Vertex>();
      for (int y = 0; y < height; y++) {
        temp.add(new Vertex(x, y));
      }
      vertices.add(temp);
    }
    Random r = new Random();
    for (int x = 0; x < width; x++) {
      for (int y = 0; y < height; y++) {
        Vertex v = vertices.get(x).get(y);
        if (x != 0) {
          v.neighborEdges.add(new Edge(v, vertices.get(x - 1).get(y), r.nextInt(1000)));
        }
        if (x != width - 1) {
          v.neighborEdges.add(new Edge(v, vertices.get(x + 1).get(y), r.nextInt(1000)));
        }
        if (y != 0) {
          v.neighborEdges.add(new Edge(v, vertices.get(x).get(y - 1), r.nextInt(1000)));
        }
        if (y != height - 1) {
          v.neighborEdges.add(new Edge(v, vertices.get(x).get(y + 1), r.nextInt(1000)));
        }
      }
    }
    return vertices;
  }
  

  // colors in each vertex given the position
  
  Color colorVertex(Vertex v) {
    if (v.x == width - 1 && v.y == height - 1) {
      return Color.magenta;
    }
    else if (v.solution) {
      return Color.blue;
    }
    else if (v.x == 0 && v.y == 0) {
      return Color.green;
    }
    else if (v.visited) {
      return Color.cyan;
    }
    else {
      return Color.white;
    }
  }
  
  //  draws out the maze with edges and vertexes
  
  public WorldScene makeScene() {
    WorldScene w = new WorldScene(width * 10, height * 10);
    for (Vertex v : vertices) {
      Color col = this.colorVertex(v);
      w.placeImageXY(new RectangleImage(10, 10, OutlineMode.SOLID, col),
          (v.x * 10) + (10 * 1 / 2), (v.y * 10) + (10 * 1 / 2));
    }
    for (Edge e : walls) {
      if (e.to.x == e.from.x) {
        w.placeImageXY(
            new RectangleImage(10, 10 / 10, OutlineMode.SOLID, Color.black),
            (e.to.x * 10) + (10 * 1 / 2),
            ((e.to.y + e.from.y) * 10 / 2) + (10 * 1 / 2));
      }
      else {
        w.placeImageXY(
            new RectangleImage(10 / 10, 10, OutlineMode.SOLID, Color.black),
            ((e.to.x + e.from.x) * 10 / 2) + (10 * 1 / 2),
            (e.to.y * 10) + (10 * 1 / 2));
      }
    }
    return w;
  }


  // gets the minimum spanning tree of a given graph with Kruskal's algorithm
  ArrayList<ArrayList<Vertex>> kruskalAlgorithm(ArrayList<ArrayList<Vertex>> graph) {
    // get all edges in the graph
    ArrayList<Edge> edgeList = createEdges(graph);
    // create an empty minimum spanning tree
    ArrayList<Edge> mst = new ArrayList<Edge>();
    // create a parent hashmap for parent finding
    HashMap<Integer, Integer> parent = new HashMap<Integer, Integer>();
    // resets all neighborEdges of each Vertex in the graph
    for (ArrayList<Vertex> row : graph) {
      for (Vertex vertex : row) {
        vertex.neighborEdges = new ArrayList<Edge>();
        parent.put(vertex.identifyVertex(), vertex.identifyVertex());
      }
    }
    // sort all edges using Collections.sort()
    Collections.sort(edgeList, new SortEdges());
    // find minimum spanning tree
    while (mst.size() < height * width - 1) {
      // get the smallest edge
      Edge edge = edgeList.get(0);
      // find the parent of the "from" vertex using union-find algorithm
      int fromParent = this.unionFind(parent, edge.from.identifyVertex());
      // find the parent of the "to" vertex using union-find algorithm
      int toParent = this.unionFind(parent, edge.to.identifyVertex());
      // if the "from" and "to" vertices are not in the same tree
      if (fromParent != toParent) {
        // add the edge to the minimum spanning tree
        mst.add(edge);
        // add the edge to the neighborEdges of the "from" vertex
        edge.from.neighborEdges.add(edge);
        // add the reverse edge to the neighborEdges of the "to" vertex
        edge.to.neighborEdges.add(new Edge(edge.to, edge.from, edge.weight));
        // merge the "from" and "to" trees by setting the "to" parent to the "from" parent
        parent.put(toParent, fromParent);
      }
      // remove the processed edge from the list of all edges
      edgeList.remove(0);
    }
    // return the original graph with updated neighborEdges
    return graph;
  }

  // used for finding values in our hashmap
  int unionFind(HashMap<Integer, Integer> hashmap, int x) {
    if (hashmap.get(x) == x) {
      return x;
    }
    else {
      return unionFind(hashmap, hashmap.get(x));
    }
  }
  

  // does tasks on every tick based on 
  // breath or depth first

  public void onTick() {
    if (bfs) {
      if (b.hasNext()) {
        b.next();
      }
    }
    if (dfs) {
      if (d.hasNext()) {
        d.next();
      }
    }
  }

  // runs a key event based on key pressed
  
  public void onKeyEvent(String key) {
    if (key.equals("b")) {
      bfs = true;
      dfs = false;
      manual = false;
      this.resetMaze();
    }
    else if (key.equals("d")) {
      dfs = true;
      bfs = false;
      manual = false;
      this.resetMaze();
    }
    else if (key.equals("m")) {
      manual = true;
      bfs = false;
      dfs = false;
      this.resetMaze();
    }
    else if (key.equals("r")) {
      this.initMaze();
    }
    // controls manual controls
    else if (manual) {
      if (m.isSolved()) {
        if (key.equals("left")) {
          m.goLeft();
        }
        else if (key.equals("up")) {
          m.goUp();
        }
        else if (key.equals("right")) {
          m.goRight();
        }
        else if (key.equals("down")) {
          m.goDown();
        }
      }
    }
  }
  

  // EXTRA CREDIT METHOD
  // resets the maze without exiting the game 

  public void resetMaze() {
    for (Vertex v : vertices) {
      v.solution = false;
      v.visited = false;
    }
    b = new BFS(vertices);
    d = new DFS(vertices);
    m = new Manual(vertices);
  }
}

// for testing purposes

class ExamplesMazeGame {

  // TESTS ALL MAZEWORLD METHODS
  
  // tests the identifyVertex() method

  void testIdentifyVertex(Tester t) {
    // Create a graph with three vertices
    Vertex v1 = new Vertex(100, 2);
    Vertex v2 = new Vertex(20, 34);
    Vertex v3 = new Vertex(2, 3);
    t.checkExpect(v1.identifyVertex(), 2100);
    t.checkExpect(v2.identifyVertex(), 34020);
    t.checkExpect(v3.identifyVertex(), 3002);
  }
  

  // test initMaze() method

  void testInitMaze(Tester t) {
    Maze m = new Maze();
    m.initMaze();
    t.checkExpect(100 * 60, m.vertices.size());
    int i = 0;
    for (int y = 0; y < 100; y++) {
      for (int x = 0; x < 60; x++) {
        Vertex v = m.vertices.get(i);
        t.checkExpect(x, v.x);
        t.checkExpect(y, v.y);
        i++;
      }
    }
  }

  // tests the makeGraph() method

  void testMakeGraph(Tester t) {
    Maze m = new Maze();
    ArrayList<ArrayList<Vertex>> grid = m.makeGraph();
    t.checkExpect(grid.size(), 100);
    t.checkExpect(grid.get(0).size(), 60);
    for (int x = 0; x < 100; x++) {
      for (int y = 0; y < 60; y++) {
        Vertex vertex = grid.get(x).get(y);
        t.checkExpect(vertex.x, x);
        t.checkExpect(vertex.y, y);
      }
    }
  }

  // tests the createWalls() method

  void testCreateWalls(Tester t) {
    Maze m = new Maze();
    ArrayList<ArrayList<Vertex>> v = new ArrayList<ArrayList<Vertex>>();
    ArrayList<Vertex> row1 = new ArrayList<Vertex>();
    ArrayList<Vertex> row2 = new ArrayList<Vertex>();
    row1.add(new Vertex(0, 0));
    row1.add(new Vertex(0, 1));
    row2.add(new Vertex(1, 0));
    row2.add(new Vertex(1, 1));
    v.add(row1);
    v.add(row2);
    ArrayList<Edge> edgeList = new ArrayList<Edge>();
    edgeList.add(new Edge(row1.get(0), row1.get(1), 1));
    edgeList.add(new Edge(row1.get(1), row2.get(1), 2));
    edgeList.add(new Edge(row2.get(1), row2.get(0), 3));
    edgeList.add(new Edge(row2.get(0), row1.get(0), 4));
    ArrayList<Edge> walls = m.createWalls(v, edgeList);
    t.checkExpect(walls.size(), 4);
    t.checkExpect((walls.get(0).from.equals(row1.get(0)) 
        && walls.get(0).to.equals(row1.get(1))), true);
    t.checkExpect((walls.get(1).from.equals(row1.get(1)) 
        && walls.get(1).to.equals(row1.get(1))), false);
    t.checkExpect((walls.get(2).from.equals(row1.get(1)) 
        && walls.get(2).to.equals(row1.get(1))), false);
  }

  // tests the createEdges() method by checking 
  // the number of edges in a graph and its properties

  void testCreateEdges(Tester t) {
    Maze m = new Maze();
    ArrayList<ArrayList<Vertex>> graph = m.makeGraph();
    ArrayList<Edge> edges = m.createEdges(graph);
    int expectedNumEdges = (90 - 1) * 60 +
        (60 - 1) * 90;
    t.checkExpect(edges.size() / 2, expectedNumEdges);
    for (Edge edge : edges) {
      int dx = Math.abs(edge.from.x - edge.to.x);
      int dy = Math.abs(edge.from.y - edge.to.y);
      t.checkExpect(((dx == 1 && dy == 0) || (dx == 0 && dy == 1)), true);
    }
  }

  // tests the unionFind() method

  void testUnionFind(Tester t) {
    Maze m = new Maze();
    // create a graph with three vertices
    Vertex v1 = new Vertex(0, 0);
    Vertex v2 = new Vertex(0, 1);
    Vertex v3 = new Vertex(1, 0);
    // create an empty parent map
    HashMap<Integer, Integer> parent = new HashMap<Integer, Integer>();
    // put each vertex in its own set
    parent.put(v1.identifyVertex(), v1.identifyVertex());
    parent.put(v2.identifyVertex(), v2.identifyVertex());
    parent.put(v3.identifyVertex(), v3.identifyVertex());
    // test that the find method returns the correct parent for each vertex
    t.checkExpect(m.unionFind(parent, v1.identifyVertex()), v1.identifyVertex());
    t.checkExpect(m.unionFind(parent, v2.identifyVertex()), v2.identifyVertex());
    t.checkExpect(m.unionFind(parent, v3.identifyVertex()), v3.identifyVertex());
  }

  // tests the kruskalAlgorithm() method

  void testKruskalAlgorithm(Tester t) {
    Maze m = new Maze();
    // create a test graph
    ArrayList<ArrayList<Vertex>> testGraph = new ArrayList<ArrayList<Vertex>>();
    ArrayList<Vertex> row1 = new ArrayList<Vertex>(Arrays.asList(new Vertex(0, 0), 
        new Vertex(1, 0), new Vertex(2, 0)));
    ArrayList<Vertex> row2 = new ArrayList<Vertex>(Arrays.asList(new Vertex(0, 1), 
        new Vertex(1, 1), new Vertex(2, 1)));
    ArrayList<Vertex> row3 = new ArrayList<Vertex>(Arrays.asList(new Vertex(0, 2), 
        new Vertex(1, 2), new Vertex(2, 2)));
    testGraph.add(row1);
    testGraph.add(row2);
    testGraph.add(row3);
    // Add some edges to the test graph
    testGraph.get(0).get(0).neighborEdges.add(new Edge(testGraph.get(0).get(0), 
        testGraph.get(0).get(1), 1));
    testGraph.get(0).get(0).neighborEdges.add(new Edge(testGraph.get(0).get(0), 
        testGraph.get(1).get(0), 2));
    testGraph.get(1).get(0).neighborEdges.add(new Edge(testGraph.get(1).get(0), 
        testGraph.get(2).get(0), 3));
    testGraph.get(1).get(0).neighborEdges.add(new Edge(testGraph.get(1).get(0), 
        testGraph.get(1).get(1), 4));
    testGraph.get(1).get(0).neighborEdges.add(new Edge(testGraph.get(1).get(0), 
        testGraph.get(0).get(0), 8));
    testGraph.get(1).get(1).neighborEdges.add(new Edge(testGraph.get(1).get(1), 
        testGraph.get(2).get(1), 6));
    testGraph.get(1).get(1).neighborEdges.add(new Edge(testGraph.get(1).get(1), 
        testGraph.get(1).get(2), 7));
    testGraph.get(1).get(1).neighborEdges.add(new Edge(testGraph.get(1).get(1), 
        testGraph.get(0).get(1), 8));
    testGraph.get(2).get(0).neighborEdges.add(new Edge(testGraph.get(2).get(0), 
        testGraph.get(2).get(1), 9));
    testGraph.get(2).get(1).neighborEdges.add(new Edge(testGraph.get(2).get(1), 
        testGraph.get(2).get(2), 10));
    // run the Kruskal's algorithm
    ArrayList<ArrayList<Vertex>> mst = m.kruskalAlgorithm(testGraph);
    // checks properties of the tree after
    t.checkExpect(mst.size(), 6);
    t.checkExpect(mst.contains(testGraph.get(0)), true);
    t.checkExpect(mst.contains(testGraph.get(1)), true);
    t.checkExpect(mst.contains(testGraph.get(2)), true);
    t.checkExpect(testGraph.get(0).get(0).neighborEdges.contains(new Edge(
        testGraph.get(0).get(0), testGraph.get(0).get(1), 1)), true);
  }

  // tests the colorVertex() method

  void testColorVertex(Tester t) {
    Maze m = new Maze();
    Vertex v1 = new Vertex(0, 0);
    Vertex v2 = new Vertex(99, 59);
    Vertex v3 = new Vertex(2, 2);

    t.checkExpect(m.colorVertex(v1), Color.green);
    t.checkExpect(m.colorVertex(v2), Color.magenta);
    t.checkExpect(m.colorVertex(v3), Color.white);
  }

  // tests the makeScene() method by making a maze
  // and comparing it with the method call

  void testMakeScene(Tester t) {
    Maze m = new Maze();
    WorldScene w = new WorldScene(500, 300);
    m.initMaze();
    ArrayList<Vertex> vertices = m.vertices;
    ArrayList<Edge> walls = m.walls;
    for (Vertex v : vertices) {
      Color colour = m.colorVertex(v);
      w.placeImageXY(new RectangleImage(8, 8, OutlineMode.SOLID, colour),
          (v.x * 8) + (8 * 1 / 2), (v.y * 8) + (8 * 1 / 2));
    }
    for (Edge e : walls) {
      if (e.to.x == e.from.x) {
        w.placeImageXY(
            new RectangleImage(8, 1, OutlineMode.SOLID, Color.black),
            (e.to.x * 8) + (8 * 1 / 2),
            ((e.to.y + e.from.y) * 8 / 2) + (8 * 1 / 2));
      }
      else {
        w.placeImageXY(
            new RectangleImage(1, 8, OutlineMode.SOLID, Color.black),
            ((e.to.x + e.from.x) * 8 / 2) + (8 * 1 / 2),
            (e.to.y * 8) + (8 * 1 / 2));
      }

    }
    t.checkExpect(w.equals(m.makeScene()), true);
  }
  
  // tests the onTick() method
  
  void testOnTick(Tester t) { 
    Maze m = new Maze();
    t.checkExpect(m.b, false);
    m.onTick();
    t.checkExpect(m.b.hasNext(), true);
    t.checkExpect(m.b, true);
    m.b.next();
    t.checkExpect(m.d, false);
    m.onTick();
    t.checkExpect(m.d.hasNext(), true);
    t.checkExpect(m.d, true);
    
  }
  
  
  // tests the onKey() method

  void testOnKeyEvent(Tester t) {
    Maze m = new Maze();
    t.checkExpect(m.b, false);
    m.onKeyEvent("b");
    t.checkExpect(m.b, true);
    m.onKeyEvent("r");
    t.checkExpect(m.d, false);
    m.onKeyEvent("d");
    t.checkExpect(m.d, true);
    m.onKeyEvent("r");
    t.checkExpect(m.m, false);
    m.onKeyEvent("m");
    t.checkExpect(m.m, true);
  }
  
  // runs the maze game by calling bigBang
  
  void testGame(Tester t) {
    Maze m = new Maze();
    m.bigBang(60 * 10, 30 * 10, 0.005);
  }
  
  // TESTS FOR BFS, DFS, and MANUAL METHODS
  
  // tests the hasNext() method
  
  void testHasNext(Tester t) {
    Maze m = new Maze();
    ArrayList<Vertex> vertices = m.vertices;
    ArrayList<Vertex> empty = new ArrayList<Vertex>();
    // breadth first
    BFS bt = new BFS(vertices);
    BFS bf = new BFS(empty);
    // depth first
    DFS dt = new DFS(vertices);
    DFS df = new DFS(empty);
    t.checkExpect(bt.hasNext(), true);
    t.checkExpect(bf.hasNext(), false);
    
    t.checkExpect(dt.hasNext(), true);
    t.checkExpect(df.hasNext(), false);
    
  }
  
  // tests the next() method for both breadth and depth first search
  
  void testNext(Tester t) {
    
    // TESTS FOR BREADTH FIRST 
    
    // create a test graph with 3 vertices and 2 edges
    Vertex v1 = new Vertex(0, 0);
    Vertex v2 = new Vertex(1, 0);
    Vertex v3 = new Vertex(0, 1);

    Edge e1 = new Edge(v1, v2, 1);
    Edge e2 = new Edge(v1, v3, 1);

    v1.neighborEdges.add(e1);
    v1.neighborEdges.add(e2);

    ArrayList<Vertex> vertices = new ArrayList<>();
    vertices.add(v1);
    vertices.add(v2);
    vertices.add(v3);

    BFS bfs = new BFS(vertices);

    // perform one step of the BFS by calling next
    bfs.next();

    // verify that the vertices are marked as part of the solution
    t.checkExpect(v2.solution, true);
    t.checkExpect(v3.solution, true);
    t.checkExpect(v1.solution, false); 

    // verify that the cameFromEdge map has been properly populated
    t.checkExpect(bfs.cameFromEdge.get(v3.identifyVertex()), v1);
    t.checkExpect(bfs.cameFromEdge.get(v2.identifyVertex()), v1);
    
    // TESTS FOR DEPTH FIRST 
    
    DFS dfs = new DFS(vertices);

    // perform one step of the BFS by calling next
    dfs.next();

    // verify that the vertices are marked as part of the solution
    t.checkExpect(v2.solution, true);
    t.checkExpect(v3.solution, true);
    t.checkExpect(v1.solution, false); 

    // verify that the cameFromEdge map has been properly populated
    t.checkExpect(dfs.cameFromEdge.get(v3.identifyVertex()), v1);
    t.checkExpect(dfs.cameFromEdge.get(v2.identifyVertex()), v1);
  }
  
  // TESTS THE METHODS FOR MANUAL CLASS
  
  // tests the isSolved() method
  
  void testIsSolved(Tester t) {
    Maze m = new Maze();
    ArrayList<Vertex> vertices = m.vertices;
    ArrayList<Vertex> empty = new ArrayList<Vertex>();
    // Manual
    Manual ma = new Manual(vertices);
    Manual ma2 = new Manual(empty);
    t.checkExpect(ma.isSolved(), true);
    t.checkExpect(ma2.isSolved(), false);
  }
  
  // tests the move() method
  
  void testMove(Tester t) {
    Vertex v1 = new Vertex(0, 0);
    Vertex v2 = new Vertex(1, 0);
    Vertex v3 = new Vertex(0, 1);

    Edge e1 = new Edge(v1, v2, 1);
    Edge e2 = new Edge(v1, v3, 1);

    v1.neighborEdges.add(e1);
    v1.neighborEdges.add(e2);

    ArrayList<Vertex> vertices = new ArrayList<>();
    vertices.add(v1);
    vertices.add(v2);
    vertices.add(v3);
    // Manual
    Manual ma = new Manual(vertices);
    t.checkExpect(ma.proceed(false, e2).equals(v1), false);
    t.checkExpect(ma.proceed(false, e2).equals(v2), true);
    t.checkExpect(ma.proceed(false, e1).equals(v1), false);  
    
  }
  
  // tests the goLeft(), goRight(), goUp(), and goDown() methods
  // in the manual class
  
  void testGoMethods(Tester t) {
    Vertex v1 = new Vertex(0, 0);
    Vertex v2 = new Vertex(1, 0);
    Vertex v3 = new Vertex(0, 1);
    Vertex v4 = new Vertex(-1, 0);
    Vertex v5 = new Vertex(0, -1);
    Edge e1 = new Edge(v1, v2, 1);
    Edge e2 = new Edge(v1, v3, 1);
    v1.neighborEdges.add(e1);
    v1.neighborEdges.add(e2);
    ArrayList<Vertex> vertices = new ArrayList<>();
    vertices.add(v1);
    vertices.add(v2);
    vertices.add(v3);
    Manual ma = new Manual(vertices);
    t.checkExpect(ma.goLeft(), v4);
    t.checkExpect(ma.goRight(), v2);
    t.checkExpect(ma.goUp(), v3);
    t.checkExpect(ma.goDown(), v4);
  }
  
  // tests the reconstruct() method
  
  void testReconstruct(Tester t) {
    
    // makes example maze
    
    Vertex v1 = new Vertex(0, 0);
    Vertex v2 = new Vertex(1, 0);
    Vertex v3 = new Vertex(0, 1);
    Vertex v4 = new Vertex(-1, 0);
    Vertex v5 = new Vertex(0, -1);
    Edge e1 = new Edge(v1, v2, 1);
    Edge e2 = new Edge(v1, v3, 1);
    v1.neighborEdges.add(e1);
    v1.neighborEdges.add(e2);
    ArrayList<Vertex> vertices = new ArrayList<>();
    vertices.add(v1);
    vertices.add(v2);
    vertices.add(v3);
    
    // creates hashmap
    HashMap<Integer, Vertex> cameFromEdge = new HashMap<>();
    cameFromEdge.put(v2.identifyVertex(), v1);
    cameFromEdge.put(v3.identifyVertex(), v2);
    cameFromEdge.put(v4.identifyVertex(), v3);

    Manual ma = new Manual(vertices);
    // calls method
    ma.reconstruct(cameFromEdge, v4);
    
    // checks that solution is set for all verticies
    t.checkExpect(v1.solution, true);
    t.checkExpect(v2.solution, true);
    t.checkExpect(v3.solution, true);
    t.checkExpect(v4.solution, true);
    
    // not called for v5;
    t.checkExpect(v5.solution, false);

  }
  
  
  
  
  // TESTS ALL METHODS FOR DEQUE, NODE, AND SENTINEL CLASSES
  // CODE/TESTS TAKEN FROM DEQUES.JAVA FROM HW 8

  // void method to initialize deques

  void initDeque() {
    deque1 = new Deque<String>();

    Sentinel<String> sentinel2 = new Sentinel<String>();
    Node<String> filler2 = new Node<String>("filler");
    Node<String> node22 = new Node<String>("bcd", filler2, filler2);
    Node<String> node21 = new Node<String>("abc", node22, sentinel2);
    Node<String> node24 = new Node<String>("def", sentinel2, filler2);
    Node<String> node23 = new Node<String>("cde", node24, node22);
    deque2 = new Deque<String>(sentinel2);

    Sentinel<Integer> sentinel3 = new Sentinel<Integer>();
    Node<Integer> filler3 = new Node<Integer>(0);
    Node<Integer> node32 = new Node<Integer>(74, filler3, filler3);
    Node<Integer> node31 = new Node<Integer>(89, node32, sentinel3);
    Node<Integer> node34 = new Node<Integer>(31, sentinel3, filler3);
    Node<Integer> node33 = new Node<Integer>(57, node34, node32);
    deque3 = new Deque<Integer>(sentinel3);
  }


  // add at head examples

  Sentinel<String> sentinel11 = new Sentinel<String>();
  Node<String> nodeAAH1 = new Node<String>("test", sentinel11, sentinel11);
  Deque<String> deque11 = new Deque<String>(sentinel11);

  Sentinel<String> sentinel21 = new Sentinel<String>();
  Node<String> fillerAAH1 = new Node<String>("filler");
  Node<String> nodeAAH0 = new Node<String>("aaa", fillerAAH1, sentinel21);
  Node<String> nodeAAH3 = new Node<String>("bcd", fillerAAH1, fillerAAH1);
  Node<String> nodeAAH2 = new Node<String>("abc", nodeAAH3, nodeAAH0);
  Node<String> nodeAAH5 = new Node<String>("def", sentinel21, fillerAAH1);
  Node<String> nodeAAH4 = new Node<String>("cde", nodeAAH5, nodeAAH3);
  Deque<String> deque21 = new Deque<String>(sentinel21);

  Sentinel<Integer> sentinel31 = new Sentinel<Integer>();
  Node<Integer> fillerAAH2 = new Node<Integer>(0);
  Node<Integer> nodeAAH6 = new Node<Integer>(16, fillerAAH2, sentinel31);
  Node<Integer> nodeAAH8 = new Node<Integer>(74, fillerAAH2, fillerAAH2);
  Node<Integer> nodeAAH7 = new Node<Integer>(89, nodeAAH8, nodeAAH6);
  Node<Integer> nodeAAH10 = new Node<Integer>(31, sentinel31, fillerAAH2);
  Node<Integer> nodeAAH9 = new Node<Integer>(57, nodeAAH10, nodeAAH8);
  Deque<Integer> deque31 = new Deque<Integer>(sentinel31);


  // add at tail examples

  Sentinel<String> sentinel12 = new Sentinel<String>();
  Node<String> nodeAAT1 = new Node<String>("test", sentinel12, sentinel12);
  Deque<String> deque12 = new Deque<String>(sentinel12);

  Sentinel<String> sentinel22 = new Sentinel<String>();
  Node<String> fillerAAT1 = new Node<String>("filler");
  Node<String> nodeAAT3 = new Node<String>("bcd", fillerAAT1, fillerAAT1);
  Node<String> nodeAAT2 = new Node<String>("abc", nodeAAT3, sentinel22);
  Node<String> nodeAAT5 = new Node<String>("def", fillerAAT1, fillerAAT1);
  Node<String> nodeAAT4 = new Node<String>("cde", nodeAAT5, nodeAAT3);
  Node<String> nodeAAT6 = new Node<String>("gih", sentinel22, nodeAAT5);
  Deque<String> deque22 = new Deque<String>(sentinel22);

  Sentinel<Integer> sentinel32 = new Sentinel<Integer>();
  Node<Integer> fillerAAT2 = new Node<Integer>(0);
  Node<Integer> nodeAAT8 = new Node<Integer>(74, fillerAAT2, fillerAAT2);
  Node<Integer> nodeAAT7 = new Node<Integer>(89, nodeAAT8, sentinel32);
  Node<Integer> nodeAAT10 = new Node<Integer>(31, fillerAAT2, fillerAAT2);
  Node<Integer> nodeAAT9 = new Node<Integer>(57, nodeAAT10, nodeAAT8);
  Node<Integer> nodeAAT11 = new Node<Integer>(16, sentinel32, nodeAAT10);
  Deque<Integer> deque32 = new Deque<Integer>(sentinel32);

  //remove from head examples

  Sentinel<String> sentinel23 = new Sentinel<String>();
  Node<String> fillerRFH = new Node<String>("filler");
  Node<String> nodeRFH2 = new Node<String>("bcd", fillerRFH, sentinel23);
  Node<String> nodeRFH4 = new Node<String>("def", sentinel23, fillerRFH);
  Node<String> nodeRFH3 = new Node<String>("cde", nodeRFH4, nodeRFH2);
  Deque<String> deque23 = new Deque<String>(sentinel23);

  Sentinel<Integer> sentinel33 = new Sentinel<Integer>();
  Node<Integer> fillerRFH2 = new Node<Integer>(0);
  Node<Integer> nodeRFH6 = new Node<Integer>(74, fillerRFH2, sentinel33);
  Node<Integer> nodeRFH8 = new Node<Integer>(31, sentinel33, fillerRFH2);
  Node<Integer> nodeRFH7 = new Node<Integer>(57, nodeRFH8, nodeRFH6);
  Deque<Integer> deque33 = new Deque<Integer>(sentinel33);

  // remove from tail examples

  Sentinel<String> sentinel24 = new Sentinel<String>();
  Node<String> fillerRFT = new Node<String>("filler");
  Node<String> nodeRFT2 = new Node<String>("bcd", fillerRFT, fillerRFT);
  Node<String> nodeRFT1 = new Node<String>("abc", nodeRFT2, sentinel24);
  Node<String> nodeRFT3 = new Node<String>("cde", sentinel24, nodeRFT2);
  Deque<String> deque24 = new Deque<String>(sentinel24);

  Sentinel<Integer> sentinel34 = new Sentinel<Integer>();
  Node<Integer> fillerRFT2 = new Node<Integer>(0);
  Node<Integer> nodeRFT5 = new Node<Integer>(74, fillerRFT2, fillerRFT2);
  Node<Integer> nodeRFT4 = new Node<Integer>(89, nodeRFT5, sentinel34);
  Node<Integer> nodeRFT6 = new Node<Integer>(57, sentinel34, nodeRFT5);
  Deque<Integer> deque34 = new Deque<Integer>(sentinel34);

  Deque<String> deque1 = new Deque<String>();
  Deque<String> deque2 = new Deque<String>();
  Deque<Integer> deque3 = new Deque<Integer>();

  // tests the size() method

  boolean testSize(Tester t) {
    this.initDeque();
    return t.checkExpect(this.deque1.size(), 0)
        && t.checkExpect(this.deque2.size(), 4)
        && t.checkExpect(this.deque3.size(), 4)
        && t.checkExpect(this.sentinel12.size(), 0);
  }

  // tests the addAtHead() method

  boolean testAddAtHead(Tester t) {
    this.initDeque();
    this.deque1.addAtHead("test");
    this.deque2.addAtHead("aaa");
    this.deque3.addAtHead(16);
    return t.checkExpect(this.deque1, this.deque11)
        && t.checkExpect(this.deque2, this.deque21)
        && t.checkExpect(this.deque3, this.deque31);
  }

  // tests the addAtTail() method

  boolean testAddAtTail(Tester t) {
    this.initDeque();
    this.deque1.addAtTail("test");
    this.deque2.addAtTail("gih");
    this.deque3.addAtTail(16);
    return t.checkExpect(this.deque1, this.deque12)
        && t.checkExpect(this.deque2, this.deque22)
        && t.checkExpect(this.deque3, this.deque32);
  }

  // tests the removeFromHead() method

  void testRemoveFromHead(Tester t) {
    this.initDeque();
    this.deque2.removeFromHead();
    this.deque3.removeFromHead();
    t.checkExpect(deque2, deque23);
    t.checkExpect(deque3, deque33);
    t.checkException(new RuntimeException(
        "Cannot remove a Sentinel"), this.sentinel34, "removeFromHead");
  }

  // tests the removeFromTail() method

  void testRemoveFromTail(Tester t) {
    this.initDeque();
    this.deque2.removeFromTail();
    this.deque3.removeFromTail();
    t.checkExpect(this.deque2, this.deque24);
    t.checkExpect(this.deque3, this.deque34);
    t.checkException(new RuntimeException(
        "Cannot remove a Sentinel"), this.sentinel34, "removeFromTail"); 
  }

  // tests the removeNode() method

  void testRemoveNode(Tester t) {
    this.initDeque();
    this.deque12.removeNode(this.nodeAAH0);
    t.checkExpect(this.deque12.header.next.next, this.nodeAAH3);
    this.deque22.removeNode(this.fillerAAT1);
    t.checkExpect(this.deque22.header, this.sentinel22);
  }

  // tests the find() method

  boolean testFind(Tester t) {
    this.initDeque();
    return t.checkExpect(this.sentinel21.find(new StringCheck("abc")), this.nodeAAH2)
        && t.checkExpect(this.sentinel22.find(new StringCheck("gih")), this.nodeAAT6)
        && t.checkExpect(this.deque12.find(new StringCheck("filler")), this.nodeAAT5);
  }
  
  // tests the getData() method
  
  void testGetData(Tester t) {
    this.initDeque();
    t.checkExpect(this.nodeAAH0.getData(), "aaa");
    t.checkExpect(this.nodeAAT1.getData(), "test");
    t.checkExpect(this.nodeRFT3.getData(), "cde");
  }
  
  // tests the hasNextVal() method
  
  void testHasNextVal(Tester t) {
    this.initDeque();
    t.checkExpect(this.fillerRFT.hasNextVal(), true);
    t.checkExpect(this.nodeAAH0.hasNextVal(), false);
    t.checkExpect(this.nodeAAH1.hasNextVal(), true);
  }
  
  // tests the isSentinel() method
  
  void testIsSentinel(Tester t) {
    this.initDeque();
    t.checkExpect(this.sentinel11.isSentinel(), true);
    t.checkExpect(this.nodeRFT5.isSentinel(), false);
    t.checkExpect(this.nodeAAT4.isSentinel(), false);
  }

  // TESTS FOR ALL METHODS IN QUEUE CLASS
  // tests isEmpty(), remove(), and add()

  Queue<String> queue = new Queue<String>();

  void testQueue(Tester t) {
    this.initDeque();
    queue.list = deque1; 
    t.checkExpect(queue.isEmpty(), false);
    t.checkExpect(queue.remove(), "harish");
    queue.add("harish");
    t.checkExpect(queue.remove(), "fundies");
    queue.empty();
    t.checkExpect(queue.isEmpty(), true);
  }
}