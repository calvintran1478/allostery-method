require "./graph"

# Maintains k trees in sorted order from greatest total weight to least.
class TreeList
  def initialize(@k : Int32, @num_nodes : Int32)
    # Determine graph size based on number of nodes
    graph_size = Graph.get_required_buffer_size(@num_nodes)

    # Allocate memory for storing tree data
    total_buffer = LibC.malloc((k + 1) * sizeof(Graph) + (k + 1) * graph_size)
    if total_buffer.null?
      puts "Error allocating memory for tree list"
      exit 1
    end

    @trees = total_buffer.as(Pointer(Graph))
    @tree_buffer = (total_buffer.as(Pointer(Graph)) + (k + 1)).as(Pointer(Float64))

    # Initialize trees to use the allocated buffers
    (k+1).times do |i|
      buffer = (@tree_buffer.as(UInt8*) + graph_size).as(Float64*)
      @trees[i] = Graph.new(@num_nodes, buffer)
      @trees[i].weight = Float64::INFINITY
    end

    @length = k
  end

  # Returns the internal buffer used for storing tree data.
  #
  # Should only be accessed after work with the tree list is done in order to
  # free its memory.
  def buffer : Void*
    @trees.as(Void*)
  end

  # Returns a reference to the tree whose contents should be written to for the
  # next add.
  #
  # This should be done by copying an existing graph into the one provided. (No
  # guarantees are made to the validity of the returned graph, so edges should
  # not be added or removed before the graph is overriden).
  #
  # This mechanism mainly serves to avoid allocating new graph objects by
  # reusing ones that are unneeded.
  def get_temp_tree : Graph
    @trees[@k]
  end

  # Adds temp tree to the list if it is not already in the list and may be one
  # of the k least-weighted trees in the collection.
  #
  # ```
  # tree_lst = TreeList.new(2)
  #
  # G1 = Graph.new(3)
  # G1.add_edge(0, 1, weight: 3)
  # G1.add_edge(1, 2, weight: 1)
  #
  # G1.copy_to(tree_lst.get_temp_tree)
  # tree_lst.add_temp_tree
  #
  # G2 = Graph.new(3)
  # G2.add_edge(0, 2, weight: 2)
  #
  # G2.copy_to(tree_lst.get_temp_tree)
  # tree_lst.add_temp_tree
  #
  # print(tree_lst.select_min == G2) # => true
  # tree_lst.remove_min
  # print(tree_lst.select_min == G1) # => true
  # ```
  def add_temp_tree : Nil
    if @trees[@k].weight <= @trees[0].weight
      # Search for placement index
      placement_index = 0
      found = false
      while placement_index < @length && @trees[@k].weight <= @trees[placement_index].weight && !found
        if @trees[@k] == @trees[placement_index]
          found = true
        end
        placement_index += 1
      end
      placement_index -= 1 if placement_index > 0

      # Update list of trees if the new tree is not already in the list
      if !found
        temp = @trees[0]
        placement_index.times do |i|
          @trees[i] = @trees[i+1]
        end
        @trees[placement_index] = @trees[@k]
        @trees[@k] = temp
      end
    end
  end

  # Removes the least-weighted tree currently held in this list.
  def remove_min : Nil
    @length -= 1
  end

  # Returns the least-weighted tree in this list.
  def select_min : Graph
    @trees[@length - 1]
  end
end
