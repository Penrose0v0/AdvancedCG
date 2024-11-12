#pragma once

// Slightly modified version of "GLSLProgramObject.h/cpp" 
// in the book "More OpenGL Game Programming."

#include <GL/glew.h>
#include <string>
#include <map>
#include <vector>

class GLSLShaderObject;

class GLSLProgramObject 
{
public:
	GLSLProgramObject();
	GLSLProgramObject(const std::string &shader, GLenum shaderType=GL_VERTEX_SHADER_ARB);
	GLSLProgramObject(const std::string &vertexShader, const std::string &fragmentShader);
	~GLSLProgramObject();
  
	void attach(GLSLShaderObject &shader);
	void attach(GLSLShaderObject *shader);
  
	void attachShaderSourceFile(const std::string &filename, GLenum shaderType=GL_VERTEX_SHADER_ARB );
	void attachShaderCodeString(const std::string &code, GLenum shaderType=GL_VERTEX_SHADER_ARB );

	void setGeometryInputOutputTypes(GLint inputType, GLint outputType);

	void detach(GLSLShaderObject &shader);
	void detach(GLSLShaderObject *shader);
	void link();
	void use() const;
	void disable() const;

	bool linkSucceeded() const { return isLinked; }

	// explicitly specify argument types, otherwise we cannot handle array types ("*v").

	void sendUniform1f(const std::string &name, float x);
	void sendUniform1fv(const std::string &name, const float * const v, int count = 1);
	void sendUniform2f(const std::string &name, float x, float y);
	void sendUniform2fv(const std::string &name, const float * const v, int count = 1);
	void sendUniform3f(const std::string &name, float x, float y, float z);
	void sendUniform3fv(const std::string &name, const float * const v, int count = 1);
	void sendUniform4f(const std::string &name, float x, float y, float z, float w);
	void sendUniform4fv(const std::string &name, const float * const v, int count = 1);
	void sendUniform1i(const std::string &name, int x);
	void sendUniformliv(const std::string &name, const int * const v, int count = 1);
	void sendUniform2i(const std::string &name, int x, int y);
	void sendUniform2iv(const std::string &name, const int * const v, int count = 1);
	void sendUniform3i(const std::string &name, int x, int y, int z);
	void sendUniform3iv(const std::string &name, const int * const v, int count = 1);
	void sendUniform4i(const std::string &name, int x, int y, int z, int w);
	void sendUniform4iv(const std::string &name, const int * const v, int count = 1);
	void sendUniformMatrix2fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix3fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix4fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix2x3fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix2x4fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix3x2fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix3x4fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix4x2fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix4x3fv(const std::string &name, const float * const m, int count = 1, bool trans = false);
	void sendUniform1ui(const std::string &name, GLuint x);
	void sendUniform1uiv(const std::string &name, GLuint *v, int count=1);
	void sendUniform2ui(const std::string &name, GLuint x, GLuint y);
	void sendUniform2uiv(const std::string &name, GLuint *v, int count=1);
	void sendUniform3ui(const std::string &name, GLuint x, GLuint y, GLuint z);
	void sendUniform3uiv(const std::string &name, GLuint *v, int count=1);
	void sendUniform4ui(const std::string &name, GLuint x, GLuint y, GLuint z, GLuint w);
	void sendUniform4uiv(const std::string &name, GLuint *v, int count=1);

	void sendUniform1f(GLuint location, float x);
	void sendUniform1fv(GLuint location, const float * const v, int count = 1);
	void sendUniform2f(GLuint location, float x, float y);
	void sendUniform2fv(GLuint location, const float * const v, int count = 1);
	void sendUniform3f(GLuint location, float x, float y, float z);
	void sendUniform3fv(GLuint location, const float * const v, int count = 1);
	void sendUniform4f(GLuint location, float x, float y, float z, float w);
	void sendUniform4fv(GLuint location, const float * const v, int count = 1);
	void sendUniform1i(GLuint location, int x);
	void sendUniformliv(GLuint location, const int * const v, int count = 1);
	void sendUniform2i(GLuint location, int x, int y);
	void sendUniform2iv(GLuint location, const int * const v, int count = 1);
	void sendUniform3i(GLuint location, int x, int y, int z);
	void sendUniform3iv(GLuint location, const int * const v, int count = 1);
	void sendUniform4i(GLuint location, int x, int y, int z, int w);
	void sendUniform4iv(GLuint location, const int * const v, int count = 1);
	void sendUniformMatrix2fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix3fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix4fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix2x3fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix2x4fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix3x2fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix3x4fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix4x2fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniformMatrix4x3fv(GLuint location, const float * const m, int count = 1, bool trans = false);
	void sendUniform1ui(GLuint location, GLuint x);
	void sendUniform1uiv(GLuint location, const GLuint * const v, int count = 1);
	void sendUniform2ui(GLuint location, GLuint x, GLuint y);
	void sendUniform2uiv(GLuint location, const GLuint * const v, int count = 1);
	void sendUniform3ui(GLuint location, GLuint x, GLuint y, GLuint z);
	void sendUniform3uiv(GLuint location, const GLuint * const v, int count = 1);
	void sendUniform4ui(GLuint location, GLuint x, GLuint y, GLuint z, GLuint w);
	void sendUniform4uiv(GLuint location, const GLuint * const v, int count = 1);

	GLuint getUniformLocation(const std::string &name) const;
	void setAttributeLocation(const std::string &name, GLuint location);
	GLuint getAttributeLocation(const std::string &name) const;

	void getProgramLog(std::string &log) const;
	void printProgramLog() const;
	GLuint getHandle() const;
	void getParameter(GLenum param, GLint *data) const;
	void validate() const;
	bool IsValidProgram() const;
	void getAttributeInfo(GLuint location, std::string &name, GLenum &type, GLint &size) const;
	void getUniformInfo(GLuint location, std::string &name, GLenum &datatype, GLint &size) const;
	void getAttachedShaders(std::vector<GLuint> &shaderhandles);
	void getAttachedShaders(std::vector<GLSLShaderObject> &shaders);

private:
	int getLoc(const std::string &name);

	//GLSLShaderObject *vertex_, *fragment_;
	std::vector<GLSLShaderObject *> shaders_;
	bool isLinked;
	GLuint handle_;
	typedef std::map<std::string,int> uniformmap;
	uniformmap uniforms_;	
};

