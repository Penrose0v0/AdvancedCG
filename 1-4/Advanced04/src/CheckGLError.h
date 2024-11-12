#pragma once

#define CHECK_GL_ERROR CheckGLError(__FUNCTION__, __FILE__, __LINE__)

void CheckGLError(const char* func, const char* file, int line);