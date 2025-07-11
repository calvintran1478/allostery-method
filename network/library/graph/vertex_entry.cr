
class Graph
  class VertexEntry
    getter vertex : Int32
    property key : Float64
    property parent : Int32

    def initialize(@vertex : Int32, @key : Float64, @parent : Int32)
    end
  end
end
