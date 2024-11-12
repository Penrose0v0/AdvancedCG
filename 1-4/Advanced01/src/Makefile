advanced01: AbstractScene.o CheckGLError.o GLSLProgramObject.o GLSLShaderObject.o Image2OGLTexture.o Scene01Checker2D.o Scene02ImageSmoothing.o Scene03WaveAnimation.o Scene04PseudoNormal.o Scene05EnvironmentMapping.o TriMesh.o arcball_camera.o imgui.o imgui_draw.o imgui_impl_glfw.o imgui_impl_opengl3.o imgui_tables.o imgui_widgets.o main.o tinyfiledialogs.o
	g++ -o advanced01 AbstractScene.o CheckGLError.o GLSLProgramObject.o GLSLShaderObject.o Image2OGLTexture.o Scene01Checker2D.o Scene02ImageSmoothing.o Scene03WaveAnimation.o Scene04PseudoNormal.o Scene05EnvironmentMapping.o TriMesh.o arcball_camera.o imgui.o imgui_draw.o imgui_impl_glfw.o imgui_impl_opengl3.o imgui_tables.o imgui_widgets.o main.o tinyfiledialogs.o -lglfw -lGLEW -lGL -lIL -lILU -lILUT
.cpp.o:
	g++ -c $< -O3 -I../../include -std=c++11
run: advanced01
	./advanced01
clean:
	rm -f *.o advanced01
