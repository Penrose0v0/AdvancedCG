#pragma once

#include <GL/glew.h>

//bool Image2OGLTexture(const char* filename, GLuint &texID, int &width, int &height, GLenum texTarget = GL_TEXTURE_2D, GLenum filterParam = GL_LINEAR, GLenum clampParam = GL_CLAMP);	// GL_CLAMP is deprecated after OpenGL 3.0
bool Image2OGLTexture(const char* filename, GLuint& texID, int& width, int& height, GLenum texTarget = GL_TEXTURE_2D, GLenum filterParam = GL_LINEAR, GLenum clampParam = GL_CLAMP_TO_EDGE);
