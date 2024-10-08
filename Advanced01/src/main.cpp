#include <iostream>
#include <cmath>
#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>

// for cross-platform file dialogs
#include "tinyfiledialogs.h"

// for DevIL, image I/O library
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

#include "Scene01Checker2D.h"
#include "Scene02ImageSmoothing.h"
#include "Scene03WaveAnimation.h"
#include "Scene04PseudoNormal.h"
#include "Scene05EnvironmentMapping.h"

using namespace std;

// scene entries

struct SceneEntry
{
	string name;
	void (*Init)();
	void (*ReloadShaders)();
	void (*Draw)();
	void (*Cursor)(GLFWwindow* window, double xpos, double ypos);
	void (*Mouse)(GLFWwindow* window, int button, int action, int mods);
	void (*Resize)(GLFWwindow* window, int w, int h);
	void (*ImGui)();
	void (*Destroy)();
};

#define SceneRegister(name, className) { name, (className::Init), (className::ReloadShaders), (className::Draw), (className::Cursor), \
										(className::Mouse), (className::Resize), (className::ImGui), (className::Destroy) }

SceneEntry g_SceneEntries[] = {
	SceneRegister("Scene 01: Checker 2D", Scene01Checker2D),
	SceneRegister("Scene 02: Image Smoothing", Scene02ImageSmoothing),
	SceneRegister("Scene 03: Wave Animation", Scene03WaveAnimation),
	SceneRegister("Scene 04: Pseudo Normal", Scene04PseudoNormal),
	SceneRegister("Scene 05: Environment Mapping", Scene05EnvironmentMapping)
};

#undef SceneRegister

const int g_NumSceneEntries = sizeof(g_SceneEntries) / sizeof(SceneEntry);
int g_SceneIndex = 0;

glm::mat4 g_ProjMatrix;

// callbacks

void errorCallback(int error, const char* description)
{
	cout << __FUNCTION__ << ": " << description << endl;
}

void cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
	g_SceneEntries[g_SceneIndex].Cursor(window, xpos, ypos);
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	g_SceneEntries[g_SceneIndex].Mouse(window, button, action, mods);
}

void keyboardCallback(GLFWwindow* window, int key, int scanCode, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case GLFW_KEY_ESCAPE:
		case GLFW_KEY_Q:
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;
		}
	}
}

void resizeCallback(GLFWwindow* window, int w, int h)
{
	if (h < 1) h = 1;
	g_SceneEntries[g_SceneIndex].Resize(window, w, h);

	const float aspect = w / float(h);
	g_ProjMatrix = glm::perspective(45.f, aspect, 0.01f, 100.f);
}

#include <glm/gtc/type_ptr.hpp>

void init()
{
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	cout << "Vendor: " << glGetString(GL_VENDOR) << endl;
	cout << "Renderer: " << glGetString(GL_RENDERER) << endl;

	if (glewInit() != GLEW_OK)
	{
		cerr << __FUNCTION__ << ": glewInit() failed" << endl;
	}

	// initialize DevIL (image I/O library)
	ilInit();
	iluInit();
	ilutInit();
	ilutRenderer(ILUT_OPENGL);

#ifdef _WIN32
	// disable unicode in file dialogs
	tinyfd_winUtf8 = 0;
#endif

	for (int i = 0; i < g_NumSceneEntries; i++)
		g_SceneEntries[i].Init();
}

int main(int argc, char **argv)
{
	if (!glfwInit()) return 1;

	glfwSetErrorCallback(errorCallback);

#ifdef __APPLE__
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);
#endif

	GLFWwindow* pWindow = glfwCreateWindow(800, 800, "Advanced CG 01: Assignments", NULL, NULL);

	glfwMakeContextCurrent(pWindow);
	glewExperimental = GL_TRUE;
	//glfwSwapInterval(0);	// try to turn off vertical sync

	init();

	glfwSetCursorPosCallback(pWindow, cursorPosCallback);
	glfwSetMouseButtonCallback(pWindow, mouseButtonCallback);
	glfwSetKeyCallback(pWindow, keyboardCallback);
	glfwSetWindowSizeCallback(pWindow, resizeCallback);

	resizeCallback(pWindow, 800, 800);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	// imgui init
	ImGui_ImplGlfw_InitForOpenGL(pWindow, true);
	ImGui_ImplOpenGL3_Init();

	glClearColor(0.45f, 0.55f, 0.60f, 1.00f);

	while (!glfwWindowShouldClose(pWindow))
	{
		glfwPollEvents();

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		g_SceneEntries[g_SceneIndex].Draw();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// gui
		{
			ImGui::Begin("ImGui Window");

			ImGui::Text("Framerate: %.3f ms/frame (%.1f fps)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

			if (ImGui::BeginListBox("Scenes", ImVec2(0, ImGui::GetTextLineHeightWithSpacing() * g_NumSceneEntries * 1.05)))	// 1.05 for padding
			{
				for (int i = 0; i < g_NumSceneEntries; i++)
				{
					const bool is_selected = (g_SceneIndex == i);
					if (ImGui::Selectable(g_SceneEntries[i].name.c_str(), is_selected))
					{
						g_SceneIndex = i;
						//cerr << "Selected: " << g_SceneEntries[i].name << endl;
					}

					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndListBox();
			}

			ImGui::Separator();

			g_SceneEntries[g_SceneIndex].ImGui();

			ImGui::End();
		}

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(pWindow);
	}

	for (int i = 0; i < g_NumSceneEntries; i++)
		g_SceneEntries[i].Destroy();

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();;

	glfwTerminate();

	return 0;
}
