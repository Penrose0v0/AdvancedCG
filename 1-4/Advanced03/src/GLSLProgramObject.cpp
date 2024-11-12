#include "GLSLProgramObject.h"
#include "GLSLShaderObject.h"
#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>

// Slightly modified version of "GLSLProgramObjectObject.h/cpp" 
// in the book "More OpenGL Game Programming."

using namespace std;

// Public part :

GLSLProgramObject::GLSLProgramObject()
  //:vertex_(NULL), fragment_(NULL)
  : isLinked(false), handle_(0)
{
	handle_ = glCreateProgram();
}

GLSLProgramObject::GLSLProgramObject(const string &shader, unsigned int shaderType)
  : isLinked(false)
{
	handle_ = glCreateProgram();

	shaders_.push_back( new GLSLShaderObject(shader, shaderType) );
	attach( shaders_.back() );

	link();
}

GLSLProgramObject::GLSLProgramObject(const string &vertexShader, const string &fragmentShader)
  : isLinked(false)
{
	handle_ = glCreateProgram();

	shaders_.push_back( new GLSLShaderObject(vertexShader, GL_VERTEX_SHADER_ARB) );
	attach( shaders_.back() );
	shaders_.push_back( new GLSLShaderObject(fragmentShader, GL_FRAGMENT_SHADER_ARB) );
	attach( shaders_.back() );

	link();
}

GLSLProgramObject::~GLSLProgramObject()
{
	for (int i=0; i<(int)shaders_.size(); i++) delete shaders_[i];
	shaders_.clear();

	glDeleteProgram(handle_);
}

GLuint GLSLProgramObject::getHandle() const
{
	return handle_;
}

void GLSLProgramObject::getProgramLog(string &log) const
{
	GLchar *debug;
	GLint debugLength;
	getParameter(GL_INFO_LOG_LENGTH, &debugLength);

	debug = new GLcharARB[debugLength];
	glGetProgramInfoLog(handle_, debugLength, &debugLength, debug);

	//cout << debug;
	log.append(debug,0,debugLength);
	delete [] debug;
}

void GLSLProgramObject::printProgramLog() const
{
	GLchar *debug;
	GLint debugLength;
	getParameter(GL_INFO_LOG_LENGTH, &debugLength);

	debug = new GLcharARB[debugLength];
	glGetProgramInfoLog(handle_, debugLength, &debugLength, debug);

	cerr << debug;

	delete [] debug;
}

void GLSLProgramObject::validate() const
{
	glValidateProgram(handle_);
}

bool GLSLProgramObject::IsValidProgram() const
{
	GLint status;
	getParameter(GL_VALIDATE_STATUS,&status);
	if(!status)
		return false;
	else
		return true;
}

void GLSLProgramObject::getParameter(GLenum param, GLint *data) const
{
	glGetProgramiv(handle_, param, data);
}

void GLSLProgramObject::attach(GLSLShaderObject &shader)
{
	glAttachShader(handle_, shader.getHandle());
}

void GLSLProgramObject::attach(GLSLShaderObject *shader)
{
	attach(*shader);
}

void GLSLProgramObject::detach(GLSLShaderObject &shader)
{
	glDetachShader(handle_,shader.getHandle());
}

void GLSLProgramObject::detach(GLSLShaderObject *shader)
{
	detach(*shader);
}

void GLSLProgramObject::attachShaderSourceFile( const string &filename, GLenum shaderType/*=GL_VERTEX_SHADER_ARB */ )
{
	shaders_.push_back( new GLSLShaderObject(filename, shaderType) );
	attach( shaders_.back() );
}

void GLSLProgramObject::attachShaderCodeString( const string &code, GLenum shaderType/*=GL_VERTEX_SHADER_ARB */ )
{
	GLSLShaderObject *pShader = new GLSLShaderObject(shaderType);
	pShader->setShaderSource( code );
	pShader->compile();

	shaders_.push_back( pShader );
	attach( shaders_.back() );
}

void GLSLProgramObject::setGeometryInputOutputTypes(GLint inputType, GLint outputType)
{
	glProgramParameteriARB(handle_, GL_GEOMETRY_INPUT_TYPE, inputType);
	glProgramParameteriARB(handle_, GL_GEOMETRY_OUTPUT_TYPE, outputType);
}

void GLSLProgramObject::link()
{
	glLinkProgram(handle_);
	uniforms_.clear();

	GLint linked;
	glGetProgramiv(handle_, GL_LINK_STATUS, &linked);

	isLinked = (linked == GL_TRUE);
}

void GLSLProgramObject::use() const
{
	if (isLinked)
		glUseProgram(handle_);
}

void GLSLProgramObject::disable() const
{
	glUseProgram(0);
}

void GLSLProgramObject::sendUniform1f(const string &name, float x) { glUniform1f(getLoc(name), x); }
void GLSLProgramObject::sendUniform1fv(const string &name, const float * const v, int count) { glUniform1fv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform2f(const string &name, float x, float y) { glUniform2f(getLoc(name), x, y); } 
void GLSLProgramObject::sendUniform2fv(const string &name, const float * const v, int count) { glUniform2fv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform3f(const string &name, float x, float y, float z) { glUniform3f(getLoc(name), x, y, z); }
void GLSLProgramObject::sendUniform3fv(const string &name, const float * const v, int count) { glUniform3fv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform4f(const string &name, float x, float y, float z, float w) { glUniform4f(getLoc(name), x, y, z, w); }
void GLSLProgramObject::sendUniform4fv(const string &name, const float * const v, int count) { glUniform4fv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform1i(const string &name, int x) { glUniform1i(getLoc(name), x); }
void GLSLProgramObject::sendUniformliv(const string &name, const int * const v, int count) { glUniform1iv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform2i(const string &name, int x, int y) { glUniform2i(getLoc(name), x, y); }
void GLSLProgramObject::sendUniform2iv(const string &name, const int * const v, int count) { glUniform2iv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform3i(const string &name, int x, int y, int z) { glUniform3i(getLoc(name), x, y, z); }
void GLSLProgramObject::sendUniform3iv(const string &name, const int * const v, int count) { glUniform3iv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform4i(const string &name, int x, int y, int z, int w) { glUniform4i(getLoc(name), x, y, z, w); }
void GLSLProgramObject::sendUniform4iv(const string &name, const int * const v, int count) { glUniform4iv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniformMatrix2fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix2fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix3fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix3fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix4fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix4fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix2x3fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix2x3fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix2x4fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix2x4fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix3x2fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix3x2fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix3x4fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix3x4fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix4x2fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix4x2fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniformMatrix4x3fv(const string &name, const float * const m, int count, bool trans) { glUniformMatrix4x3fv(getLoc(name), count, trans, m); }
void GLSLProgramObject::sendUniform1ui(const string &name, GLuint x) { glUniform1ui(getLoc(name), x); }
void GLSLProgramObject::sendUniform1uiv(const string &name, GLuint *v, int count) { glUniform1uiv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform2ui(const string &name, GLuint x, GLuint y) { glUniform2ui(getLoc(name), x, y); }
void GLSLProgramObject::sendUniform2uiv(const string &name, GLuint *v, int count) { glUniform2uiv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform3ui(const string &name, GLuint x, GLuint y, GLuint z) { glUniform3ui(getLoc(name), x, y, z); }
void GLSLProgramObject::sendUniform3uiv(const string &name, GLuint *v, int count) { glUniform3uiv(getLoc(name), count, v); }
void GLSLProgramObject::sendUniform4ui(const string &name, GLuint x, GLuint y, GLuint z, GLuint w) { glUniform4ui(getLoc(name), x, y, z, w); }
void GLSLProgramObject::sendUniform4uiv(const string &name, GLuint *v, int count) { glUniform4uiv(getLoc(name), count, v); }

void GLSLProgramObject::sendUniform1f(GLuint location, float x) { glUniform1f(location, x); }
void GLSLProgramObject::sendUniform1fv(GLuint location, const float * const v, int count) { glUniform1fv(location, count, v); }
void GLSLProgramObject::sendUniform2f(GLuint location, float x, float y) { glUniform2f(location, x, y); } 
void GLSLProgramObject::sendUniform2fv(GLuint location, const float * const v, int count) { glUniform2fv(location, count, v); }
void GLSLProgramObject::sendUniform3f(GLuint location, float x, float y, float z) { glUniform3f(location, x, y, z); }
void GLSLProgramObject::sendUniform3fv(GLuint location, const float * const v, int count) { glUniform3fv(location, count, v); }
void GLSLProgramObject::sendUniform4f(GLuint location, float x, float y, float z, float w) { glUniform4f(location, x, y, z, w); }
void GLSLProgramObject::sendUniform4fv(GLuint location, const float * const v, int count) { glUniform4fv(location, count, v); }
void GLSLProgramObject::sendUniform1i(GLuint location, int x) { glUniform1i(location, x); }
void GLSLProgramObject::sendUniformliv(GLuint location, const int * const v, int count) { glUniform1iv(location, count, v); }
void GLSLProgramObject::sendUniform2i(GLuint location, int x, int y) { glUniform2i(location, x, y); }
void GLSLProgramObject::sendUniform2iv(GLuint location, const int * const v, int count) { glUniform2iv(location, count, v); }
void GLSLProgramObject::sendUniform3i(GLuint location, int x, int y, int z) { glUniform3i(location, x, y, z); }
void GLSLProgramObject::sendUniform3iv(GLuint location, const int * const v, int count) { glUniform3iv(location, count, v); }
void GLSLProgramObject::sendUniform4i(GLuint location, int x, int y, int z, int w) { glUniform4i(location, x, y, z, w); }
void GLSLProgramObject::sendUniform4iv(GLuint location, const int * const v, int count) { glUniform4iv(location, count, v); }
void GLSLProgramObject::sendUniformMatrix2fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix2fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix3fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix3fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix4fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix4fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix2x3fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix2x3fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix2x4fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix2x4fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix3x2fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix3x2fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix3x4fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix3x4fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix4x2fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix4x2fv(location, count, trans, m); }
void GLSLProgramObject::sendUniformMatrix4x3fv(GLuint location, const float * const m, int count, bool trans) { glUniformMatrix4x3fv(location, count, trans, m); }
void GLSLProgramObject::sendUniform1ui(GLuint location, GLuint x) { glUniform1ui(location, x); }
void GLSLProgramObject::sendUniform1uiv(GLuint location, const GLuint * const v, int count) { glUniform1uiv(location, count, v); }
void GLSLProgramObject::sendUniform2ui(GLuint location, GLuint x, GLuint y) { glUniform2ui(location, x, y); }
void GLSLProgramObject::sendUniform2uiv(GLuint location, const GLuint * const v, int count) { glUniform2uiv(location, count, v); }
void GLSLProgramObject::sendUniform3ui(GLuint location, GLuint x, GLuint y, GLuint z) { glUniform3ui(location, x, y, z); }
void GLSLProgramObject::sendUniform3uiv(GLuint location, const GLuint * const v, int count) { glUniform3uiv(location, count, v); }
void GLSLProgramObject::sendUniform4ui(GLuint location, GLuint x, GLuint y, GLuint z, GLuint w) { glUniform4ui(location, x, y, z, w); }
void GLSLProgramObject::sendUniform4uiv(GLuint location, const GLuint * const v, int count) { glUniform4uiv(location, count, v); }

GLuint GLSLProgramObject::getUniformLocation(const string &name) const
{
	return glGetUniformLocation(handle_,name.c_str());
}

void GLSLProgramObject::setAttributeLocation(const string &name, GLuint location)
{
	glBindAttribLocation(handle_,location, name.c_str());
}

GLuint GLSLProgramObject::getAttributeLocation(const string &name) const
{
	return glGetAttribLocation(handle_,name.c_str());	
}

void GLSLProgramObject::getAttributeInfo(GLuint location, string &name, GLenum &datatype, GLint &size) const
{
	// first we need to query the length of the longest name
	GLsizei length;
	getParameter(GL_ACTIVE_ATTRIBUTE_MAX_LENGTH,(GLint*)&length);

	vector<char> nametmp(length);
	GLsizei copiedsize;
	glGetActiveAttrib(handle_,location,length,&copiedsize,&size,&datatype,&nametmp[0]);

	if(0 == copiedsize)
		return;

	name.assign(&nametmp[0]);
}

void GLSLProgramObject::getUniformInfo(GLuint location, string &name, GLenum &datatype, GLint &size) const
{
	// first we need to query the length of the longest name
	GLint length;
	getParameter(GL_ACTIVE_UNIFORM_MAX_LENGTH,&length);

	vector<char> nametmp(length);
	GLsizei copiedsize;
	glGetActiveUniform(handle_,location,length,&copiedsize,&size,&datatype,&nametmp[0]);

	if(0 == copiedsize)
		return;

	name.assign(&nametmp[0]);
}

void GLSLProgramObject::getAttachedShaders(vector<GLuint> &shaderhandles)
{
	GLint number;
	getParameter(GL_ATTACHED_SHADERS,&number);

	shaderhandles.resize(number);
	glGetAttachedShaders(handle_,number,NULL, &shaderhandles[0]);
}

namespace {

	struct shaderconstuctor
	{
		vector<GLSLShaderObject> &shaders_;
		shaderconstuctor(vector<GLSLShaderObject> &shaders) : shaders_(shaders)
		{
		}
		void operator()(vector<GLuint>::value_type data)
		{
			shaders_.push_back(GLSLShaderObject(data));
		}
	};
}

void GLSLProgramObject::getAttachedShaders(vector<GLSLShaderObject> &shaders)
{
	vector<GLuint> shaderhandles;
	getAttachedShaders(shaderhandles);

	for_each(shaderhandles.begin(),shaderhandles.end(),shaderconstuctor(shaders));
}

// Private part :

int GLSLProgramObject::getLoc(const string &name)
{
	int loc = -1;
	uniformmap::const_iterator it = uniforms_.find(name);
	if(it == uniforms_.end())
	{
		loc = glGetUniformLocation(handle_, name.c_str());
		if(loc == -1)
		{
			//string s;
			//s.append(name);
			//s.append(" - is not a valid uniform variable name");
			//throw logic_error(s);
			cerr << "Error: \"" << name << "\" is not a valid uniform variable name" << endl;
			return -1;
 		}
		uniforms_.insert(uniformmap::value_type(name,loc));
	}
	else
		loc = it->second;

	return loc;
}
