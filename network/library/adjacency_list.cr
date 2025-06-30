
# An adjacency list representing edge connections within a graph.
#
# The main use for this data structure is to quickly determine the neighbours of
# each node when the graph is sparse (e.g., when the graph is a tree).
struct AdjacencyList

  # Initializes an adjacency list representing a graph with the given number of
  # nodes.
  #
  # The buffer is expected to have a capacity of num_nodes ** 2. The buffer can
  # contain any content prior to calling this method. (This content is cleared
  # upon initialization.)
  def initialize(@num_nodes : Int32, @buffer : Pointer(Int32))
    clear
  end

  # Adds edge to the adjacency list.
  #
  # Precondition: Edge (i, j) does not exist in the list.
  def add_edge(i : Int32, j : Int32) : Nil
    # Update neighbours of vertex i
    row_i = @buffer + (i * @num_nodes)
    row_i[row_i[0] + 1] = j
    row_i[0] += 1

    # Update neighbours of vertex j
    row_j = @buffer + (j * @num_nodes)
    num_neighbours =
    row_j[row_j[0] + 1] = i
    row_j[0] += 1
  end

  # Removes edge from the adjacency list.
  #
  # Precondition: Edge (i, j) exists in the list.
  def remove_edge(i : Int32, j : Int32) : Nil
    # Remove j from vertex i's neighbours
    row_i = @buffer + (i * @num_nodes)
    swap_index = -1
    curr_index = 1
    while swap_index == -1
      swap_index = curr_index if row_i[curr_index] == j
      curr_index += 1
    end

    row_i[swap_index] = row_i[row_i[0]]
    row_i[0] -= 1

    # Remove i from vertex j's neighbours
    row_j = @buffer + (j * @num_nodes)
    swap_index = -1
    curr_index = 1
    while swap_index == -1
      swap_index = curr_index if row_j[curr_index] == i
      curr_index += 1
    end

    row_j[swap_index] = row_j[row_j[0]]
    row_j[0] -= 1
  end

  @[AlwaysInline]
  def [](index : Int32) : Slice(Int32)
    Slice.new(@buffer + (index * @num_nodes + 1), @buffer[index * @num_nodes])
  end

  # Iterates over each of the given node's neighbours
  def each(node : Int32, &block : Int32 ->) : Nil
    num_neighbours = @buffer[node * @num_nodes]
    neighbours = @buffer + (node * @num_nodes + 1)

    num_neighbours.times do |i|
      yield neighbours[i]
    end
  end

  # Clears all edges in the list
  def clear : Nil
    @num_nodes.times { |i| @buffer[i * @num_nodes] = 0 }
  end
end
