#include <GL/glew.h>
#include <iostream>
#include "CheckGLError.h"

// see https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glGetError.xhtml

void CheckGLError(const char *func, const char *file, int line)
{
	GLenum errCode = glGetError();

	if (errCode != GL_NO_ERROR)
	{
		std::cerr << func << ": OpenGL Error: ";

		switch (errCode)
		{
		case GL_INVALID_ENUM:
			std::cerr << "GL_INVALID_ENUM: An unacceptable value is specified for an enumerated argument";
			break;
		case GL_INVALID_VALUE:
			std::cerr << "GL_INVALID_VALUE: A numeric argument is out of range";
			break;
		case GL_INVALID_OPERATION:
			std::cerr << "GL_INVALID_OPERATION: The specified operation is not allowed in the current state";
			break;
		case GL_STACK_OVERFLOW:
			std::cerr << "GL_STACK_OVERFLOW: This command would cause a stack overflow";
			break;
		case GL_STACK_UNDERFLOW:
			std::cerr << "GL_STACK_UNDERFLOW: This command would cause a stack underflow";
			break;
		case GL_OUT_OF_MEMORY:
			std::cerr << "GL_OUT_OF_MEMORY: There is not enough memory left to execute the command";
			break;
		case GL_TABLE_TOO_LARGE:
			std::cerr << "GL_TABLE_TOO_LARGE: The specified table exceeds the implementation's maximum supported table size";
			break;
		default:
			std::cerr << "Unknown command, error code = " << std::showbase << std::hex << errCode;
			break;
		}

		std::cerr << " (file: " << file << " at line " << line << std::endl;
	}
}
