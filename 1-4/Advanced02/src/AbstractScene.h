#pragma once

#include "GLSLProgramObject.h"
#include "GLFW/glfw3.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/ext.hpp"	// for glm::to_string()

class AbstractScene
{
public:
	static void Resize(GLFWwindow* window, int w, int h)
	{
		if (h == 0) h = 1;
		glViewport(0, 0, w, h);
		s_WindowWidth = w;
		s_WindowHeight = h;
	}

protected:
	static int s_WindowWidth, s_WindowHeight;

};