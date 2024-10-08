#include "Image2OGLTexture.h"
#include <iostream>
#undef _UNICODE
#include <IL/il.h>

using namespace std;

bool Image2OGLTexture(const char* filename, GLuint& texID, int& width, int& height, GLenum texTarget, GLenum filterParam, GLenum clampParam)
{
	ILuint imgName;
	ilGenImages(1, &imgName);
	ilBindImage(imgName);

	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_LOWER_LEFT);

	if (!ilLoadImage(filename))
	{
		cerr << __FUNCTION__ << ": Error: cannot open " << filename << endl;
		ilDeleteImages(1, &imgName);
		return false;
	}

	width = ilGetInteger(IL_IMAGE_WIDTH);
	height = ilGetInteger(IL_IMAGE_HEIGHT);

	if (!texID) glGenTextures(1, &texID);
	glBindTexture(texTarget, texID);
	glTexParameteri(texTarget, GL_TEXTURE_MAG_FILTER, filterParam);
	glTexParameteri(texTarget, GL_TEXTURE_MIN_FILTER, filterParam);
	glTexParameteri(texTarget, GL_TEXTURE_WRAP_S, clampParam);
	glTexParameteri(texTarget, GL_TEXTURE_WRAP_T, clampParam);
	glTexImage2D(texTarget, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, ilGetData());
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);	// deprecated from OpenGL 3.0
	glBindTexture(texTarget, 0);

	ilDeleteImages(1, &imgName);

	cerr << __FUNCTION__ << ": file loaded: " << filename << " (" << width << "x" << height << ")" << endl;

	return true;
}
