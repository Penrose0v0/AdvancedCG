#include "GLSLShaderObject.h"
#include <string>
#include <fstream>
#include <iostream>
#include <string.h>	// for Ubuntu (memset)

// Slightly modified version of "GLSLShader.h/cpp" 
// in the book "More OpenGL Game Programming."

using namespace std;

// Public part :

// too fat constructor ... 
GLSLShaderObject::GLSLShaderObject(const std::string &filename, GLenum shaderType) : compiled_(false)
{
	shaderType_ = shaderType;
	handle_ = glCreateShader(shaderType);

	const char * source = readShader(filename);
	
	if ( 0 == source)
	{
		cerr << __FUNCTION__ << ": shader not constructed" << endl;
	}
	else
	{
		glShaderSource(handle_, 1, static_cast<const GLchar**>(&source), NULL);

		compile();

		if ( ! isCompiled())
		{
			cerr << __FUNCTION__ << ": cannot compile shader: " << filename << endl;
			string log;
			getShaderLog(log);
			cerr << log;
		}

		delete [] source;
	}
}

GLSLShaderObject::GLSLShaderObject(GLenum shaderType)
{
	shaderType_ = shaderType;
	handle_ = glCreateShader(shaderType);
}


GLSLShaderObject::~GLSLShaderObject()
{
	glDeleteShader(handle_);
}

void GLSLShaderObject::compile()
{
	GLint compiled;

	glCompileShader(handle_);
	getParameter(GL_COMPILE_STATUS, &compiled);
	compiled_ = (compiled == GL_TRUE);
}

void GLSLShaderObject::getShaderLog(std::string &log) const
{
	GLchar *debug;
	GLint debugLength;
	getParameter(GL_INFO_LOG_LENGTH, &debugLength);
	  
	debug = new GLchar[debugLength];
	glGetShaderInfoLog(handle_, debugLength, &debugLength, debug);
		
	//cout << debug;
	log.append(debug,0,debugLength);
	delete [] debug;
}

void GLSLShaderObject::getShaderSource(std::string &shader) const
{
	GLchar *source;
	GLint sourcelen;
	getParameter(GL_SHADER_SOURCE_LENGTH, &sourcelen);

	source = new GLchar[sourcelen];
	glGetShaderSource(handle_,sourcelen, &sourcelen, source);
	shader.append(source,0,sourcelen);
	delete [] source;
}

GLuint GLSLShaderObject::getHandle() const
{
	return handle_;
}

void GLSLShaderObject::getParameter(GLenum param, GLint *data) const
{
	glGetShaderiv(handle_, param, data);
}

void GLSLShaderObject::setShaderSource(const std::string &code)
{
	filename_ = "(source string)";

	const char *source = code.c_str();
	glShaderSource(handle_,1,static_cast<const GLchar**>(&source),NULL);
}

// Private part :

char *GLSLShaderObject::readShader(const std::string &filename)
{
	filename_ = filename;

	ifstream temp(filename.c_str() );
	if ( ! temp)
	{
		cerr << __FUNCTION__ << ": cannot open file: " << filename << endl;
		return 0;
	}

	int count = 0;
	char *buf;

	temp.seekg(0, ios::end);
	count = temp.tellg();

	buf = new char[count + 1];
	memset(buf,0,count);
	temp.seekg(0, ios::beg);
	temp.read(buf, count);
	buf[count] = 0;
	temp.close();

	return buf;
}
