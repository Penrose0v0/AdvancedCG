TARGET=advanced04

$(TARGET): BlinnPhongRenderer.o CatmullClarkSubdivision.o CheckGLError.o EnvironmentMap.o GLSLProgramObject.o GLSLShaderObject.o HalfEdgeDataStructure.o LoopSubdivision.o PolygonMesh.o ReflectionLineRenderer.o arcball_camera.o imgui.o imgui_demo.o imgui_draw.o imgui_impl_glfw.o imgui_impl_opengl2.o imgui_tables.o imgui_widgets.o main.o tinyfiledialogs.o
	g++ -o $(TARGET) BlinnPhongRenderer.o CatmullClarkSubdivision.o CheckGLError.o EnvironmentMap.o GLSLProgramObject.o GLSLShaderObject.o HalfEdgeDataStructure.o LoopSubdivision.o PolygonMesh.o ReflectionLineRenderer.o arcball_camera.o imgui.o imgui_demo.o imgui_draw.o imgui_impl_glfw.o imgui_impl_opengl2.o imgui_tables.o imgui_widgets.o main.o tinyfiledialogs.o -lglfw -lGLEW -lGL -lIL -lILU -lILUT
.cpp.o:
	g++ -c $< -O3 -I../../include -std=c++11
run: $(TARGET)
	./$(TARGET)
clean:
	rm -f *.o $(TARGET)
