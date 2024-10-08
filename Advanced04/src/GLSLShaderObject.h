#pragma once

// Slightly modified version of "GLSLShader.h/cpp" 
// in the book "More OpenGL Game Programming."

#include <GL/glew.h>
#include <string>

class GLSLShaderObject
{
public:
	GLSLShaderObject(const std::string &filename, GLenum shaderType = GL_VERTEX_SHADER);
	GLSLShaderObject(GLenum shaderType = GL_VERTEX_SHADER );
	~GLSLShaderObject();
	void compile();
	inline bool isCompiled() const { return compiled_; } 
	void getShaderLog(std::string &log) const;
	void getShaderSource(std::string &shader) const;
	void setShaderSource(const std::string &code);

	GLuint getHandle() const;
	void getParameter(GLenum param, GLint *data) const;

	inline std::string getFilename() const { return filename_; }
	inline GLenum getShaderType() const { return shaderType_; }

private:
	std::string filename_;	// added for debugging
	GLenum shaderType_;		// added for debugging
	char *readShader(const std::string &filename);
	bool compiled_;
	GLint handle_;
};
