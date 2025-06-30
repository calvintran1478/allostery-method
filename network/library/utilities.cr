require "csv"
require "./graph"

# Loads an undirected graph from a csv of edge weights
def read_data(filename : String) : Tuple(Graph, String)
  # Read csv file
  file_data = File.read(filename)
  csv_data = CSV.parse(file_data)
  num_nodes = csv_data.size

  # Create initial graph
  graph = Graph.new(num_nodes)

  # Check that the data represents a square matrix, symmetrical matrix
  if csv_data.any? { |row| row.size != num_nodes }
    return graph, "Error: Input CSV must represent a square matrix"
  end

  num_nodes.times do |i|
    j = i + 1
    while j < num_nodes
      edge_weight = csv_data[i][j].to_f64?
      if csv_data[i][j] != csv_data[j][i]
        return graph, "Error: Input CSV must represent a symmetric matrix"
      end
      j += 1
    end
  end

  # Add edges to the graph
  graph = Graph.new(num_nodes)
  csv_data.each_with_index do |row, i|
    j = i + 1
    while j < num_nodes
      # Check that each entry in the matrix is a positive number
      edge_weight = row[j].to_f64?
      if edge_weight.nil? || edge_weight <= 0
        return graph, "Error: CSV entries must be positive numbers"
      end

      graph.add_edge(i, j, row[j].to_f64)
      j += 1
    end
  end

  return graph, ""
end

# Saves the graph as a csv of its edge weights
def save_graph(filename : String, graph : Graph) : Nil
  File.open(filename, "w") do |file|
    graph.num_nodes.times do |i|
      edge_value = (i != 0) ? graph.get_edge(i, 0) : 0
      file << edge_value
      j = 1
      while j < graph.num_nodes
        edge_value = (i != j) ? graph.get_edge(i, j) : 0
        file << ',' << edge_value
        j += 1
      end
      file << '\n'
    end
  end
end
