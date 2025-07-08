require "./adjacency_list"
require "./tree_list"

# Type aliases
alias Edge = Tuple(Int32, Int32, Float64)
alias VertexEntry = NamedTuple(vertex: Int32, key: Float64, parent: Int32)

# Undirected graph containing a fixed number of nodes.
#
# Only accepts nonnegative edge weights.
class Graph
  getter num_nodes : Int32
  getter num_edges : Int32
  property weight : Float64
  protected getter buffer : Float64*

  # Initializes the graph with the given number of nodes.
  #
  # Accepts an optional buffer for storing edge weights and other information.
  # The user should call Graph.get_required_buffer_size to get the required
  # number of bytes for the buffer.
  #
  # ```
  # num_nodes = 10
  # buffer_size = Graph.get_required_buffer_size(num_nodes)
  # buffer = LibC.malloc(buffer_size).as(Float64*)
  # graph = Graph.new(num_nodes, buffer)
  #
  # LibC.free(graph.buffer)
  # ```
  def initialize(@num_nodes : Int32, buffer : Float64* | Nil = nil) : Nil
    # Allocate buffer for storing graph data or set buffer to the one provided
    @num_edges = (num_nodes * (num_nodes - 1)) >> 1
    @self_allocated = buffer.nil?
    @buffer = !buffer.nil? ? buffer : LibC.malloc(Graph.get_required_buffer_size(num_nodes)).as(Float64*)

    # Check if memory allocation is successful
    if @buffer.null?
      puts "Error allocating memory for graph"
      exit 1
    end

    # Calculate offsets for quick edge weight lookup
    @offset_buffer = (@buffer + @num_edges).as(Int32*)
    @offset_buffer[0] = 0
    (1...num_nodes - 1).each { |i| @offset_buffer[i] = @offset_buffer[i-1] + (num_nodes - i) }

    # Clear edge buffer and set initial weight to 0
    @num_edges.times { |i| @buffer[i] = -1 }
    @weight = 0
  end

  # Returns the number of bytes needed for constructing a graph with the given
  # number of nodes.
  #
  # Intended purpose is for manually allocating memory for a graph to avoid
  # placing added pressure on the garbage collector.
  def Graph.get_required_buffer_size(num_nodes : Int32) : Int32
    num_rows = num_nodes - 1
    num_edges = (num_nodes * (num_nodes - 1)) >> 1 # Equivalent to (num_nodes * (num_nodes - 1) / 2).to_i; optimized using bit shift operator

    (num_edges * sizeof(Float64)) + (num_rows * sizeof(Int32))
  end

  # Returns the index into self.buffer which gives the edge weight between nodes
  # i and j.
  @[AlwaysInline]
  private def get_edge_index(i : Int32, j : Int32) : Int32
    (i < j) ? @offset_buffer[i] + (j - i - 1) : @offset_buffer[j] + (i - j - 1)
  end

  # Frees internal resources used by the graph when the graph is eventually
  # garbage collected.
  def finalize
    LibC.free(@buffer) if @self_allocated
  end

  # Returns the edge weight between vertices i and j if it exists, and -1
  # otherwise.
  @[AlwaysInline]
  def get_edge(i : Int32, j : Int32) : Float64
    @buffer[get_edge_index(i, j)]
  end

  # Adds an edge between nodes i and j with the given weight.
  #
  # Assumes there is no prior edge between nodes i and j and that the weight is
  # nonnegative.
  def add_edge(i : Int32, j : Int32, weight : Float64) : Nil
    @buffer[get_edge_index(i, j)] = weight
    @weight += weight
  end

  # Removes the edge between vertices i and j if it exists.
  #
  # If no edge exists between vertices i and j no changes are made to the graph.
  def remove_edge(i : Int32, j : Int32) : Nil
    index = get_edge_index(i, j)

    @weight -= @buffer[index]
    @buffer[index] = -1
  end

  # Copies self into the given graph.
  #
  # Assumes the other graph has same number of nodes as self.
  def copy_to(other : Graph) : Nil
    @buffer.copy_to(other.buffer, @num_edges)
    other.weight = @weight
  end

  # Writes the minimum spanning tree (MST) of self to the given graph.
  #
  # Assumes the given graph has the same number of nodes as self and contains
  # no edges. MST construction is done using Prim's algorithm.
  #
  # Implementation Reference: Cormen, T.H., Leiserson, C.E., Rivest, R.L., &
  # Stein, C. (2009). Introduction to Algorithms, third edition.
  #
  # ```
  # G1 = Graph.new(3)
  # G1.add_edge(0, 1, weight: 1)
  # G1.add_edge(1, 2, weight: 2)
  # G1.add_edge(0, 2, weight: 3)
  #
  # G2 = Graph.new(3)
  # G1.write_minimum_spanning_tree(G2)
  # print(G2) # => [(0, 1, weight=1), (1, 2, weight=2)]
  # ```
  def write_minimum_spanning_tree(empty_graph : Graph) : Nil
    # Initialize MST
    mst = empty_graph
    mst.weight = 0

    # Initialize min-priority queue for quickly finding which node is closest to
    # the current tree. Also stores which connections the MST should have
    vertex_buffer = LibC.malloc(@num_nodes * sizeof(VertexEntry)).as(Pointer(VertexEntry))
    if vertex_buffer.null?
      puts "Error allocating memory for MST calculation"
      exit 1
    end

    @num_nodes.times { |i| vertex_buffer[i] = {vertex: i, key: Float64::INFINITY, parent: -1} }

    # Initialize root as starting node for MST construction
    vertex_buffer[0] = {vertex: 0, key: 0.0, parent: -1}

    # Construct MST by iteratively adding nodes which are closest to the current
    # tree
    (@num_nodes - 1).times do |i|
      # Extract min (closest node) from the priority queue
      u = vertex_buffer[i]

      # Update priority queue by finding which node is now closest
      j = i + 1
      closest_vertex_index = j
      closest_weight = vertex_buffer[j][:key]
      while j < @num_nodes
        v = vertex_buffer[j]
        if 0 < get_edge(u[:vertex], v[:vertex]) < v[:key]
          vertex_buffer[j] = {vertex: v[:vertex], key: get_edge(u[:vertex], v[:vertex]), parent: u[:vertex]}
        end

        v = vertex_buffer[j]
        if v[:key] < closest_weight
          closest_vertex = v[:vertex]
          closest_weight = v[:key]
        end

        j += 1
      end

      # Move new closest node to the front of the queue
      vertex_buffer[closest_vertex_index], vertex_buffer[i + 1] = vertex_buffer[i + 1], vertex_buffer[closest_vertex_index]
    end

    # Build MST from the determined edges
    i = 1
    while i < @num_nodes
      mst.add_edge(vertex_buffer[i][:vertex], vertex_buffer[i][:parent], get_edge(vertex_buffer[i][:vertex], vertex_buffer[i][:parent]))
      i += 1
    end

    # Free allocated memory (we don't need the vertex buffer anymore)
    LibC.free(vertex_buffer)
  end

  # Iterates over each edge in the graph and yields the two endpoints of the
  # edge along with its weight.
  #
  # ```
  # G1 = Graph.new(3)
  # G1.add_edge(1, 2, weight: 1)
  # G1.add_edge(0, 2, weight: 2)
  #
  # G1.each_edge do |i, j, weight|
  #   print({i, j, weight}) # => prints {1, 2, 1}, and then {0, 2, 2}
  # end
  # ```
  def each_edge(&block : Edge ->) : Nil
    index = 0
    (@num_nodes - 1).times do |i|
      (i + 1...@num_nodes).each do |j|
        yield ({i, j, @buffer[index]}) if @buffer[index] != -1
        index += 1
      end
    end
  end

  # Returns an adjacency list representing edge connections within this graph.
  #
  # The adjacency list buffer should be large enough to hold (N ** 2) values,
  # where N is the number of nodes in this graph.
  #
  # ```
  # G1 = Graph.new(3)
  # G1.add_edge(0, 1, weight: 1)
  # G1.add_edge(1, 2, weight: 2)
  #
  # adjacency_list_buffer = LibC.malloc((3 ** 2) * sizeof(Int32)).as(Int32*)
  # adjacency_list = G1.get_adjacency_list(adjacency_list_buffer)
  # adjacency_list.each(1) do |neighbour|
  #   print(neighbour) # => prints 0 and then 2
  # end
  #
  # LibC.free(adjacency_list_buffer)
  # ```
  def get_adjacency_list(adjacency_list_buffer : Pointer(Int32)) : AdjacencyList
    adjacency_list = AdjacencyList.new(@num_nodes, adjacency_list_buffer)
    each_edge do |i, j, weight|
      adjacency_list.add_edge(i, j)
    end

    adjacency_list
  end

  # Iterates over edges in self that are not in the provided graph and yields
  # the two endpoints of the edge along with its weight.
  #
  # Assumes the two graphs have the same number of nodes.
  #
  # ```
  # G1 = Graph.new(3)
  # G1.add_edge(0, 1, weight: 1)
  # G1.add_edge(1, 2, weight: 2)
  # G1.add_edge(0, 2, weight: 3)
  #
  # G2 = Graph.new(3)
  # G2.add_edge(1, 2, weight: 2)
  #
  # G1.get_disjoint_edges(G2) do |i, j, weight|
  #   print({i, j, weight}) # => prints {0, 1, 1}, and then {0, 2, 3}
  # end
  # ```
  def get_disjoint_edges(other : Graph, &block : Edge ->) : Nil
    each_edge do |i, j, weight|
      yield ({i, j, weight}) if other.get_edge(i, j) == -1
    end
  end

  # Finds all edges that participate in a cycle within the graph and writes them
  # to the given edge buffer.
  #
  # Assumes that there is only one cycle in the graph, the starting node
  # provided is part of this cycle, and the adjacency list accurately reflects
  # this graph's edges.
  #
  # The edge buffer should be able to hold N edges, where N is the number of
  # nodes in this graph. The frontier and parent buffers should also have a
  # capacity of N. (TODO: Needs fixing, and finish docs)
  #
  # ```
  # G1 = Graph.new(4)
  # G1.add_edge(0, 1, weight: 1)
  # G1.add_edge(1, 2, weight: 2)
  # G1.add_edge(0, 2, weight: 3)
  # G1.add_edge(2, 3, weight: 4)
  #
  #
  # ```
  def find_cycle(starting_node : Int32, adjacency_list : AdjacencyList, edge_buffer : Pointer(Edge), frontier : Pointer(Int32), parent : Pointer(Int32)) : Int32
    parent[starting_node] = -1
    frontier[0] = starting_node
    insert_pos = 1

    while parent[starting_node] == -1
      # Extract node from frontier
      curr_node = frontier[insert_pos - 1]
      insert_pos -= 1

      adjacency_list.each(curr_node) do |neighbour|
        if neighbour != parent[curr_node]
          parent[neighbour] = curr_node
          frontier[insert_pos] = neighbour
          insert_pos += 1
        end
      end
    end

    p "start writing"
    # Write solution to edge buffer
    edge_buffer[0] = {starting_node, parent[starting_node], get_edge(starting_node, parent[starting_node])}
    curr_node = parent[starting_node]
    i = 1
    while curr_node != starting_node
      if curr_node < parent[curr_node]
        edge_buffer[i] = {curr_node, parent[curr_node], get_edge(curr_node, parent[curr_node])}
      else
        edge_buffer[i] = {parent[curr_node], curr_node, get_edge(curr_node, parent[curr_node])}
      end

      curr_node = parent[curr_node]
      i += 1
    end

    i
  end

  # Returns the k least-weighted spanning trees from this graph.
  #
  # Returns a tuple of pointers that should be freed once work with trees is
  # finished.
  #
  # Implementation Reference: Amal, P. M., & KS, A. K. (2016). An algorithm for
  # kth minimum spanning tree. Electronic Notes in Discrete Mathematics, 53,
  # 343-354.
  #
  # TODO: Fix and finish docs
  #
  # ```
  #
  # ```
  def get_k_minimum_spanning_trees(k : Int32) : Tuple(Graph*, Void*)
    # Compute first MST
    tree_lst = TreeList.new(k, @num_nodes)
    write_minimum_spanning_tree(tree_lst.get_temp_tree)
    tree_lst.add_temp_tree

    k_MSTs = LibC.malloc(k * sizeof(Graph)).as(Pointer(Graph))
    if k_MSTs.null?
      puts "Error allocating memory for k least-weighted spanning trees calculation"
      exit 1
    end
    k_MSTs[0] = tree_lst.select_min

    # Allocate memory for working with edges
    buffer = LibC.malloc(@num_nodes * sizeof(Edge) + ((@num_nodes + 2) * @num_nodes) * sizeof(Int32))
    if buffer.null?
      puts "Error allocating memory for k least-weighted spanning trees calculation"
      exit 1
    end

    cycle_edges = buffer.as(Pointer(Edge))
    adjacency_list_buffer = (cycle_edges + @num_nodes).as(Pointer(Int32))
    frontier = (adjacency_list_buffer + (@num_nodes * @num_nodes))
    parent_buffer = (frontier + @num_nodes)

    i = 1
    while i < k
      tree_lst.remove_min
      adjacency_list = k_MSTs[i-1].get_adjacency_list(adjacency_list_buffer)

      get_disjoint_edges(k_MSTs[i-1]) do |test_edge|
        k_MSTs[i-1].add_edge(test_edge[0], test_edge[1], test_edge[2])
        adjacency_list.add_edge(test_edge[0], test_edge[1])

        # Find cycle in the graph formed by adding the new edge
        num_cycle_edges = k_MSTs[i-1].find_cycle(test_edge[0], adjacency_list, cycle_edges, frontier, parent_buffer)

        # Determine weight criterion for exchange edges
        weight_criterion = -1
        if i == 1
          num_cycle_edges.times do |l|
            if weight_criterion < cycle_edges[l][2] <= test_edge[2] && cycle_edges[l] != test_edge
              weight_criterion = cycle_edges[l][2]
            end
          end
        elsif
          num_cycle_edges.times do |l|
            if weight_criterion < cycle_edges[l][2] < test_edge[2] && cycle_edges[l] != test_edge
              weight_criterion = cycle_edges[l][2]
            end
          end
        end

        # Add candidate trees to tree list
        num_cycle_edges.times do |l|
          if cycle_edges[l][2] == weight_criterion
            k_MSTs[i-1].copy_to(tree_lst.get_temp_tree)
            tree_lst.get_temp_tree.remove_edge(cycle_edges[l][0], cycle_edges[l][1])
            tree_lst.add_temp_tree
          end
        end

        adjacency_list.remove_edge(test_edge[0], test_edge[1])
        k_MSTs[i-1].remove_edge(test_edge[0], test_edge[1])
      end

      # Record discovered tree
      k_MSTs[i] = tree_lst.select_min
      i += 1
    end

    # Free allocated memory
    LibC.free(buffer)

    {k_MSTs, tree_lst.buffer}
  end

  # Constructs a network by unioning the k least-weighted spanning trees of self
  # together into a single graph.
  def construct_network(k : Int32) : Graph
    # Determine k least-weighted spanning trees
    k_MSTs, residual_buffer = get_k_minimum_spanning_trees(k)

    # Combine these trees into a single network by unioning them together
    network = Graph.new(@num_nodes)
    (@num_nodes - 1).times do |i|
      (i + 1...@num_nodes).each do |j|
        edge_value = -1.0
        k.times do |l|
          edge_value = Math.max(edge_value, k_MSTs[l].get_edge(i, j))
        end
        network.add_edge(i, j, edge_value)
      end
    end

    # Free used memory
    LibC.free(residual_buffer)
    LibC.free(k_MSTs)

    network
  end

  # Returns whether the graphs have the same edges.
  #
  # Assumes that they have the same number of nodes.
  #
  # ```
  # G1 = Graph.new(3)
  # G2 = Graph.new(3)
  # G1.add_edge(0, 1, weight: 1)
  # G2.add_edge(0, 1, weight: 1)
  #
  # print(G1 == G2) # => true
  #
  # G2.add_edge(1, 2, weight: 2)
  #
  # print(G1 == G2) # => false
  # ```
  def ==(other : Graph) : Bool
    # Check that both graphs have the same weight
    return false if @weight != other.weight

    # Check that all edges are the same
    return @buffer.memcmp(other.buffer, @num_edges) == 0
  end

  # Prints edges of the graph to the given io
  #
  # ```
  # G1 = Graph.new(3)
  # G1.add_edge(0, 1, weight: 1)
  # G1.add_edge(1, 2, weight: 2)
  #
  # print(G1) # => [(0, 1, weight=1.0), (1, 2, weight=2.0)]
  # ```
  def to_s(io : IO) : Nil
    index = 0
    first = true

    io << "["
    each_edge do |i, j, weight|
      if first
        io << "(" << i << ", " << j << ", weight=" << weight << ")"
        first = false
      else
        io << ", (" << i << ", " << j << ", weight=" << weight << ")"
      end
    end
    io << "]"
  end
end
