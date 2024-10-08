#include <iostream>
#include <cmath>
#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>

#include "PathFinder.h"

// for cross-platform file dialogs
#include "tinyfiledialogs.h"

// for DevIL, image I/O library
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

#include "TriMesh.h"
#include "arcball_camera.h"

#include "DirectionalLightManager.h"
#include "MaterialParameters.h"

#include "Scene01ShadingExamples.h"
#include "Scene02ShadowMapping.h"
#include "Scene03MultipleRenderTarget.h"

using namespace std;

// main contents of the scenes

int g_WindowWidth = 800, g_WindowHeight = 800;

TriMesh g_TriMesh, g_Plane;
MaterialParameters g_Material;

glm::vec2 g_PrevMouse;
ArcballCamera g_Camera(glm::vec3(0.f, 1.f, 5.f), glm::vec3(0.f, 1.f, 0.f), glm::vec3(0, 1, 0));
glm::mat4 g_ProjMatrix = glm::perspective(45.f, g_WindowWidth / float(g_WindowHeight), 0.01f, 10.f);

bool g_DraggingMenu = false;

// scene entries

struct SceneEntry
{
	string name;
	void (*Init)();
	void (*ReloadShaders)();
	void (*Draw)();
	void (*Resize)(GLFWwindow* window, int w, int h);
	void (*ImGui)();
	void (*Destroy)();
};

#define SceneRegister(name, className) { name, (className::Init), (className::ReloadShaders), (className::Draw), \
										(className::Resize), (className::ImGui), (className::Destroy) }

SceneEntry g_SceneEntries[] = {
	SceneRegister("Scene 01: Shading Examples", Scene01ShadingExamples),
	SceneRegister("Scene 02: Shadow Mapping", Scene02ShadowMapping),
	SceneRegister("Scene 03: Multiple Render Targets", Scene03MultipleRenderTarget)
};

#undef SceneRegister

const int g_NumSceneEntries = sizeof(g_SceneEntries) / sizeof(SceneEntry);
int g_SceneIndex = 0;

// callbacks

void errorCallback(int error, const char* description)
{
	cout << __FUNCTION__ << ": " << description << endl;
}

void cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
	if (g_DraggingMenu)
		return;

	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE)
	{
		return;
	}

	glm::vec2 currPos(xpos / g_WindowWidth, (g_WindowHeight - ypos) / g_WindowHeight);

	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
	{
		g_Camera.rotate(g_PrevMouse, currPos);
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS)
	{
		g_Camera.zoom(currPos.y - g_PrevMouse.y);
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
	{
		g_Camera.pan(currPos - g_PrevMouse);
	}

	g_PrevMouse = currPos;
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		g_PrevMouse = glm::vec2(xpos / g_WindowWidth, (g_WindowHeight - ypos) / g_WindowHeight);
	}
}

void keyboardCallback(GLFWwindow* window, int key, int scanCode, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case GLFW_KEY_C:
			cerr << "Camera: " << glm::to_string(g_Camera.transform()) << endl;
			break;
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

	g_WindowWidth = w;
	g_WindowHeight = h;
	g_ProjMatrix = glm::perspective(45.f, w / float(h), 0.01f, 100.f);
}

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

	// define a plane as two triangles
	{
		const float halfSize = 100.f;
		//const float planeY = -1.f;
		const float planeY = -0.75f;
		vector<glm::vec3> planeVertices = {
			glm::vec3(-halfSize, planeY, -halfSize), glm::vec3(-halfSize, planeY, halfSize), glm::vec3(halfSize, planeY, halfSize), glm::vec3(halfSize, planeY, -halfSize)
		};
		vector<TriangleIndices> triangleIndices = { TriangleIndices(0, 1, 2), TriangleIndices(0, 2, 3) };

		g_Plane.setVertices(planeVertices);
		g_Plane.setTriangles(triangleIndices);
		g_Plane.calcVertexNormals();
		g_Plane.bakeVBOs();
	}

	// load a mesh
	{
		PathFinder finder;
		finder.addSearchPath("Resources");
		finder.addSearchPath("../Resources");
		finder.addSearchPath("../../Resources");

		g_TriMesh.loadObj(finder.find("duck.obj").c_str());
		//g_TriMesh.loadTexture(finder.find("duckCM.jpg").c_str());
	}

	DirectionalLightManager::Init();

	//glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
	//glClampColor(GL_CLAMP_VERTEX_COLOR, GL_FALSE);
	//glClampColor(GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);

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

	GLFWwindow* pWindow = glfwCreateWindow(g_WindowWidth, g_WindowHeight, "Advanced CG 02: Assignments", NULL, NULL);

	glfwMakeContextCurrent(pWindow);
	glewExperimental = GL_TRUE;
	glfwSwapInterval(0);	// try to turn off vertical sync

	init();

	glfwSetCursorPosCallback(pWindow, cursorPosCallback);
	glfwSetMouseButtonCallback(pWindow, mouseButtonCallback);
	glfwSetKeyCallback(pWindow, keyboardCallback);
	glfwSetWindowSizeCallback(pWindow, resizeCallback);

	resizeCallback(pWindow, g_WindowWidth, g_WindowHeight);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	// imgui init
	ImGui_ImplGlfw_InitForOpenGL(pWindow, true);
	ImGui_ImplOpenGL3_Init();

	glEnable(GL_DEPTH_TEST);

	//bool displayXyzAxes = false;

	while (!glfwWindowShouldClose(pWindow))
	{
		glfwPollEvents();

		glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//glMatrixMode(GL_PROJECTION);
		//glLoadMatrixf(glm::value_ptr(g_ProjMatrix));
		//glMatrixMode(GL_MODELVIEW);
		//glLoadMatrixf(glm::value_ptr(g_Camera.transform()));
		//
		//if (displayXyzAxes)
		//{
		//	glBegin(GL_LINES);
		//	glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(1, 0, 0);
		//	glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 1, 0);
		//	glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 1);
		//	glEnd();
		//}

		g_SceneEntries[g_SceneIndex].Draw();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// GUI
		{
			ImGui::Begin("ImGui Window");

			g_DraggingMenu = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);

			ImGui::Text("Framerate: %.3f ms/frame (%.1f fps)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

			if (ImGui::Button("Load Mesh"))
			{
				char const* lFilterPatterns[1] = { "*.obj" };
				const char* lTheOpenFileName = tinyfd_openFileDialog("Opening an OBJ mesh", "", 1, lFilterPatterns, "OBJ (*.obj)", 0);
				if (lTheOpenFileName)
					g_TriMesh.loadObj(lTheOpenFileName);
			}

			//ImGui::Checkbox("Display XYZ Axes", &displayXyzAxes);

			if (ImGui::BeginListBox("Scenes", ImVec2(0, ImGui::GetTextLineHeightWithSpacing() * g_NumSceneEntries * 1.05)))	// 1.05 for padding
			{
				for (int i = 0; i < g_NumSceneEntries; i++)
				{
					const bool is_selected = (g_SceneIndex == i);
					if (ImGui::Selectable(g_SceneEntries[i].name.c_str(), is_selected))
					{
						g_SceneIndex = i;
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
