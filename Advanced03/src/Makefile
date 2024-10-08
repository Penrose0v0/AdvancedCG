TARGET=advanced03

$(TARGET): CheckGLError.o EnvironmentMap.o GLSLProgramObject.o GLSLShaderObject.o GeometricObject.o Material.o PathTracer.o Scene.o Sphere.o Texture.o Triangle.o TriangleMesh.o arcball_camera.o imgui.o imgui_demo.o imgui_draw.o imgui_impl_glfw.o imgui_impl_opengl2.o imgui_tables.o imgui_widgets.o main.o tinyfiledialogs.o
	g++ -o $(TARGET) CheckGLError.o EnvironmentMap.o GLSLProgramObject.o GLSLShaderObject.o GeometricObject.o Material.o PathTracer.o Scene.o Sphere.o Texture.o Triangle.o TriangleMesh.o arcball_camera.o imgui.o imgui_demo.o imgui_draw.o imgui_impl_glfw.o imgui_impl_opengl2.o imgui_tables.o imgui_widgets.o main.o tinyfiledialogs.o -lglfw -lGLEW -lGL -lIL -lILU -lILUT -fopenmp
.cpp.o:
	g++ -c $< -O3 -I../../include -std=c++11 -fopenmp
run: $(TARGET)
	./$(TARGET)
clean:
	rm -f *.o $(TARGET)
