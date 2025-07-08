require "./graph"

# Maintains k trees in sorted order from greatest total weight to least.
#
# Internally stores a list of tuples where the first element is the tree and the
# second is its total edge sum.
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
      buffer = (@tree_buffer.as(Int8*) + graph_size).as(Float64*)
      @trees[i] = Graph.new(@num_nodes, buffer)
      @trees[i].weight = Float64::INFINITY
    end

    # Store current number of held trees (will change later)
    @length = k
  end

  # Returns the internal buffer used for storing tree data.
  #
  # Should only be accessed after work with the tree list is done in order to
  # free this memory.
  def buffer : Void*
    @trees.as(Void*)
  end

  # Returns a reference to the tree whose contents should be written to for the
  # next add
  def get_temp_tree : Graph
    @trees[@k]
  end

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
  #
  # Also reduces the capacity of the list by 1.
  def remove_min : Nil
    @length -= 1
  end

  # Returns the least-weighted tree in this list
  def select_min : Graph
    @trees[@length - 1]
  end
end
