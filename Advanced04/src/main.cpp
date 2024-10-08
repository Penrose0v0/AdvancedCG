#include <cstdio>
#include <iostream>
#include <cmath>
#include <GL/glew.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
//#include "imgui_impl_opengl3.h"
#include "GLFW/glfw3.h"
#include "PathFinder.h"
#include "CheckGLError.h"

#include "tinyfiledialogs.h"

#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

#include "PolygonMesh.h"
#include "arcball_camera.h"

#include "LoopSubdivision.h"
#include "CatmullClarkSubdivision.h"

#include "BlinnPhongRenderer.h"
#include "ReflectionLineRenderer.h"

#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "glfw3dll.lib")

#pragma comment(lib, "ILU.lib")
#pragma comment(lib, "ILUT.lib")
#pragma comment(lib, "DevIL.lib")

using namespace std;
using namespace glm;

int g_WindowWidth = 800, g_WindowHeight = 800;

GLFWwindow* g_pWindow = 0;

glm::vec2 g_PrevMouse;
ArcballCamera g_Camera(glm::vec3(1.5f), glm::vec3(0.f), glm::vec3(0, 1, 0));

glm::mat4 g_ProjMatrix = glm::perspective(45.f, g_WindowWidth / float(g_WindowHeight), 0.1f, 200.f);

PolygonMesh g_Mesh;
string g_MeshFilename = "cube.obj";

bool g_DraggingMenu = false;

// subdivision schemes

struct SubdivisionEntry
{
	string name;
	AbstractSubdivision *(*Create)();
};

#define RegisterSubsidivionScheme(name, className) { name, className::Create }

SubdivisionEntry g_SubdivisionSchemes[] = {
	RegisterSubsidivionScheme("Loop", LoopSubdivision),
	RegisterSubsidivionScheme("Catmull-Clark", CatmullClarkSubdivision)
};

#undef RegisterSubsidivionScheme

const int g_NumSubdivisionEntries = sizeof(g_SubdivisionSchemes) / sizeof(SubdivisionEntry);
int g_SubdivisionIndex = 0;

// renderers

struct RendererEntry
{
	string name;
	void (*Init)();
	AbstractRenderer* (*Create)();
	void (*ImGui)();
};

#define RegisterRenderer(name, className) { name, className::Init, className::Create, className::ImGui }

RendererEntry g_RendererEntries[] = {
	RegisterRenderer("Blinn-Phong", BlinnPhongRenderer),
	RegisterRenderer("Reflection Line", ReflectionLineRenderer)
};

#undef RegisterRenderer

const int g_NumRendererEntries = sizeof(g_RendererEntries) / sizeof(RendererEntry);
int g_RendererIndex = 0;

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

	glViewport(0, 0, w, h);
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

	ilInit();
	iluInit();
	ilutInit();
	ilutRenderer(ILUT_OPENGL);

#ifdef _WIN32
	// disable unicode in file dialogs
	tinyfd_winUtf8 = 0;
#endif

	{
		PathFinder finder;
		finder.addSearchPath("Resources");
		finder.addSearchPath("../Resources");
		finder.addSearchPath("../../Resources");

		g_MeshFilename = finder.find(g_MeshFilename);
		g_Mesh.loadObj(g_MeshFilename.c_str());
	}

	for (int ri = 0; ri < g_NumRendererEntries; ++ri)
		g_RendererEntries[ri].Init();

	glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
	glClampColor(GL_CLAMP_VERTEX_COLOR, GL_FALSE);
	glClampColor(GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);
}

int main() {

	if (!glfwInit()) return 1;

	glfwSetErrorCallback(errorCallback);

#ifdef __APPLE__
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);
#endif

	g_pWindow = glfwCreateWindow(g_WindowWidth, g_WindowHeight, "Advanced CG 04: Assignments", NULL, NULL);

	glfwMakeContextCurrent(g_pWindow);
	glewExperimental = GL_TRUE;
	//glfwSwapInterval(0);	// try to turn off vertical sync

	init();

	glfwSetCursorPosCallback(g_pWindow, cursorPosCallback);
	glfwSetMouseButtonCallback(g_pWindow, mouseButtonCallback);
	glfwSetKeyCallback(g_pWindow, keyboardCallback);
	glfwSetWindowSizeCallback(g_pWindow, resizeCallback);

	resizeCallback(g_pWindow, g_WindowWidth, g_WindowHeight);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();

	// imgui init
	ImGui_ImplGlfw_InitForOpenGL(g_pWindow, true);
	ImGui_ImplOpenGL2_Init();
	//ImGui_ImplOpenGL3_Init();

	glEnable(GL_DEPTH_TEST);

	bool displayXYZAxes = false;
	bool displayWireframe = false;

	auto pRenderer = g_RendererEntries[g_RendererIndex].Create();

	while (!glfwWindowShouldClose(g_pWindow))
	{
		glfwPollEvents();

		CHECK_GL_ERROR;

		glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_PROJECTION);
		glLoadMatrixf(glm::value_ptr(g_ProjMatrix));
		glMatrixMode(GL_MODELVIEW);
		glLoadMatrixf(glm::value_ptr(g_Camera.transform()));

		if (displayWireframe)
		{
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1, 1);
		}

		pRenderer->render(g_Mesh);

		if (displayWireframe)
		{
			glDisable(GL_POLYGON_OFFSET_FILL);

			glColor3f(0.2f, 0.2f, 0.2f);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			g_Mesh.renderMeshWithoutNormals();
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		if (displayXYZAxes)
		{
			const float axisLength = 3.f;
			glBegin(GL_LINES);
			glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(axisLength, 0, 0);
			glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, axisLength, 0);
			glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, axisLength);
			glEnd();
		}

		CHECK_GL_ERROR;

		ImGui_ImplOpenGL2_NewFrame();
		//ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// gui
		{
			ImGui::Begin("ImGui Window");

			g_DraggingMenu = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);

			ImGui::Checkbox("Display XYZ Axes", &displayXYZAxes);
			ImGui::Checkbox("Display Wireframe", &displayWireframe);

			ImGui::Separator();

			if (ImGui::Button("Load Mesh"))
			{
				char const* lFilterPatterns[1] = { "*.obj" };
				const char* lTheOpenFileName = tinyfd_openFileDialog("Opening an OBJ mesh", "", 1, lFilterPatterns, "OBJ (*.obj)", 0);
				if (lTheOpenFileName)
				{
					g_Mesh.loadObj(lTheOpenFileName);
					g_MeshFilename = string(lTheOpenFileName);
				}
			}

			ImGui::Text("# faces = %d, # verts = %d", g_Mesh.getNumFaces(), g_Mesh.getNumVertices());

			if (ImGui::Button("Reload Mesh") && g_MeshFilename != "")
			{
				g_Mesh.loadObj(g_MeshFilename.c_str());
			}

			ImGui::Separator();

			if (ImGui::BeginListBox("Subdivision", ImVec2(0, ImGui::GetTextLineHeightWithSpacing() * g_NumSubdivisionEntries * 1.05)))	// 1.05 for padding
			{
				for (int i = 0; i < g_NumSubdivisionEntries; i++)
				{
					const bool is_selected = (g_SubdivisionIndex == i);
					if (ImGui::Selectable(g_SubdivisionSchemes[i].name.c_str(), is_selected))
					{
						g_SubdivisionIndex = i;
					}

					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndListBox();
			}

			static int nSubdiv = 1;

			ImGui::SliderInt("# Subdivisions", &nSubdiv, 0, 5);

			if (ImGui::Button("Apply Subdivision"))
			{
				auto pSubdiv = g_SubdivisionSchemes[g_SubdivisionIndex].Create();
				pSubdiv->subdivide(g_Mesh, nSubdiv);
				delete pSubdiv;
			}

			ImGui::Separator();

			if (ImGui::BeginListBox("Renderers", ImVec2(0, ImGui::GetTextLineHeightWithSpacing() * g_NumRendererEntries * 1.05)))	// 1.05 for padding
			{
				for (int i = 0; i < g_NumRendererEntries; i++)
				{
					const bool is_selected = (g_RendererIndex == i);
					if (ImGui::Selectable(g_RendererEntries[i].name.c_str(), is_selected))
					{
						g_RendererIndex = i;
						delete pRenderer;
						pRenderer = g_RendererEntries[g_RendererIndex].Create();
					}

					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndListBox();
			}

			g_RendererEntries[g_RendererIndex].ImGui();

			ImGui::End();
		}

		ImGui::Render();
		ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
		//ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(g_pWindow);
	}

	delete pRenderer;

	ImGui_ImplOpenGL2_Shutdown();
	//ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();;

	glfwTerminate();

	return 0;
}