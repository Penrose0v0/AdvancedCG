//#include "stdio.h"
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

#include "arcball_camera.h"

#include "PathTracer.h"

#include "Scene.h"

#include "Sphere.h"
#include "TriangleMesh.h"

#include "PseudoNormalColorMaterial.h"
#include "DiffuseMaterial.h"
#include "BlinnPhongMaterial.h"
#include "PerfectSpecularMaterial.h"
#include "SpecularRefractionMaterial.h"

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

PathTracer g_PathTracer;

Scene g_Scene;

glm::vec2 g_PrevMouse;
//ArcballCamera g_Camera(glm::vec3(0.5f, 4.f, 5.f), glm::vec3(0.f, 1.5f, 0.f), glm::vec3(0, 1, 0));
ArcballCamera g_Camera(glm::vec3(-6, 2, 0), glm::vec3(0, 1.5, 0), glm::vec3(0, 1, 0));

glm::mat4 g_ProjMatrix = glm::perspective(45.f, g_WindowWidth / float(g_WindowHeight), 0.1f, 200.f);

bool g_KeepTracing = false;
bool g_DraggingMenu = false;
bool g_DisplayPathTracedResult = false;

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

	//g_SceneEntries[g_SceneIndex].Resize(window, w, h);

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

	glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
	glClampColor(GL_CLAMP_VERTEX_COLOR, GL_FALSE);
	glClampColor(GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);

	{
		{
			PathFinder finder;
			finder.addSearchPath("Resources");
			finder.addSearchPath("../Resources");
			finder.addSearchPath("../../Resources");

			g_Scene.loadEnvironmentMap(finder.find("sunset_fairway_2k.hdr").c_str());
		}

		// floor
		{
			DiffuseMaterial* m = DiffuseMaterial::CreateMaterial();
			m->setDiffuseCoeff(0.5f, 0.5f, 0.5f);

			const float halfSize = 3.f;

			TriangleMesh* o = TriangleMesh::CreateGeometricObject();
			o->addTriangle(Triangle(vec3(-halfSize, 0.f, halfSize), vec3(halfSize, 0.f, halfSize), vec3(halfSize, 0.f, -halfSize), m));
			o->addTriangle(Triangle(vec3(-halfSize, 0.f, halfSize), vec3(halfSize, 0.f, -halfSize), vec3(-halfSize, 0.f, -halfSize), m));
			o->bakeVBO();
			o->setMaterial(m);

			g_Scene.addObject(o);
		}

		// pyramid top
		{
			PerfectSpecularMaterial* m = PerfectSpecularMaterial::CreateMaterial();
			m->setSpecularCoeff(0.5f, 0.5f, 0.5f);
			g_Scene.addObject(Sphere::CreateGeometricObject(glm::vec3(0.f, sqrtf(2.f) + 1.f, 0.f), 1.f, m));
		}

		// front left
		{
			//DiffuseMaterial* m = DiffuseMaterial::CreateMaterial();
			BlinnPhongMaterial* m = BlinnPhongMaterial::CreateMaterial();
			m->setDiffuseCoeff(0.2f, 0.2f, 0.2f);
			m->setSpecularCoeff(0.8f, 0.2f, 0.2f);
			m->setShininess(64.f);
			g_Scene.addObject(Sphere::CreateGeometricObject(glm::vec3(-1.f, 1.f, 1.f), 1.f, m));
		}

		// front right
		{
			//DiffuseMaterial* m = DiffuseMaterial::CreateMaterial();
			//m->setDiffuseCoeff(0.5f, 1.f, 0.5f);
			SpecularRefractionMaterial* m = SpecularRefractionMaterial::CreateMaterial();
			m->setRefractionIndex(1.5f);
			m->setSpecularCoeff(0.5f, 1.f, 0.5f);
			g_Scene.addObject(Sphere::CreateGeometricObject(glm::vec3(1.f, 1.f, 1.f), 1.f, m));
		}

		// rear right
		{
			DiffuseMaterial* m = DiffuseMaterial::CreateMaterial();
			m->setDiffuseCoeff(0.5f, 1.f, 0.5f);
			g_Scene.addObject(Sphere::CreateGeometricObject(glm::vec3(-1.f, 1.f, -1.f), 1.f, m));
		}

		// rear right
		{
			DiffuseMaterial* m = DiffuseMaterial::CreateMaterial();
			m->setDiffuseCoeff(0.5f, 0.5f, 1.f);
			g_Scene.addObject(Sphere::CreateGeometricObject(glm::vec3(1.f, 1.f, -1.f), 1.f, m));
		}
	}
}

int main() {

	if (!glfwInit()) return 1;

	glfwSetErrorCallback(errorCallback);

#ifdef __APPLE__
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);
#endif

	g_pWindow = glfwCreateWindow(g_WindowWidth, g_WindowHeight, "Advanced CG 03: Assignments", NULL, NULL);

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

	while (!glfwWindowShouldClose(g_pWindow)) {
		glfwPollEvents();

		CheckGLError(__FUNCTION__, __FILE__, __LINE__);

		glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (g_DisplayPathTracedResult)
		{
			if (g_KeepTracing)
				g_PathTracer.renderScene();
			g_PathTracer.renderFrame();
		}
		else
		{
			glMatrixMode(GL_PROJECTION);
			glLoadMatrixf(glm::value_ptr(g_ProjMatrix));
			glMatrixMode(GL_MODELVIEW);
			glLoadMatrixf(glm::value_ptr(g_Camera.transform()));

			if (displayXYZAxes)
			{
				const float axisLength = 3.f;
				glBegin(GL_LINES);
				glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(axisLength, 0, 0);
				glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, axisLength, 0);
				glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, axisLength);
				glEnd();
			}

			g_Scene.drawGL();
		}

		CheckGLError(__FUNCTION__, __FILE__, __LINE__);

		ImGui_ImplOpenGL2_NewFrame();
		//ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		// gui
		{
			ImGui::Begin("ImGui Window");

			g_DraggingMenu = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);

			ImGui::Checkbox("Display XYZ Axes", &displayXYZAxes);
			//ImGui::Checkbox("Keep Tracing", &g_KeepTracing);

			if (ImGui::Button("Load Environment Map"))
			{
				char const* lFilterPatterns[] = { "*.jpg", "*.png", "*.hdr", "*.exr" };
				const char* lTheOpenFileName = tinyfd_openFileDialog(
					"Opening an environment map", "", 4, lFilterPatterns, "Environment map (*.jpg;*.png;*.hdr;*.exr)", 0);

				if (lTheOpenFileName)
					g_Scene.loadEnvironmentMap(lTheOpenFileName);
			}

			ImGui::SliderInt("Max Recursion Depth", &PathTracer::s_MaxRecursionDepth, 0, 64);
			ImGui::SliderInt("Min Recursion Depth", &PathTracer::s_MinRecursionDepth, 0, 64);
			ImGui::SliderInt("# Samples Per Pixel", &PathTracer::s_NumSamplesPerPixel, 1, 4096);
			ImGui::SliderInt("# Samples Per Update", &PathTracer::s_NumSamplesPerUpdate, 1, 4096);

			if (ImGui::Button("Render"))
			{
				g_DisplayPathTracedResult = true;
				g_PathTracer.renderScene();
			}

			ImGui::Checkbox("Display Path Traced Result", &g_DisplayPathTracedResult);

			ImGui::SliderFloat("Gamma Correction", &PathTracer::s_Gamma, 0.001f, 5.f);

			ImGui::End();
		}

		ImGui::Render();
		ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
		//ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(g_pWindow);
	}

	GeometricObject::ClearGeometricObjectCache();
	Material::ClearMaterialCache();

	ImGui_ImplOpenGL2_Shutdown();
	//ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();;

	glfwTerminate();

	return 0;
}