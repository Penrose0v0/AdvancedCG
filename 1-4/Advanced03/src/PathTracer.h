#pragma once

#include "Ray.h"
#include "ImageRect.h"
#include "glm/glm.hpp"
#include "GLSLProgramObject.h"
#include <random>

class PathTracer
{
public:
	typedef ImageRect<glm::vec3> ImageRGBf;

	static float s_Gamma;
	static int s_MaxRecursionDepth;
	static int s_MinRecursionDepth;
	static int s_NumSamplesPerPixel;
	static int s_NumSamplesPerUpdate;

	PathTracer() : m_FrameBufferTexID(0), m_RandSrc(12345), m_RandDist(0.f, 1.f), m_pGammaShader(0) {}
	~PathTracer()
	{
		if (m_FrameBufferTexID) glDeleteTextures(1, &m_FrameBufferTexID);
		if (m_pGammaShader) delete m_pGammaShader;
	}

	void renderScene();
	void renderFrame();

private:
	ImageRGBf m_FrameBuffer;
	GLuint m_FrameBufferTexID;

	std::mt19937 m_RandSrc;
	std::uniform_real_distribution<float> m_RandDist;

	GLSLProgramObject* m_pGammaShader;
	bool m_isNVIDIADriver;

	void initShader();

	glm::vec3 traceRec(const Ray& ray, int recursionDepth);

	void updateFrameBufferTexture();
	void renderIntermediateFrame();

	void calcLocalCoordinateSystem(const glm::vec3& normal, const glm::vec3& inDir, glm::vec3& xLocal, glm::vec3& yLocal, glm::vec3& zLocal) const;

	inline float frand() { return m_RandDist(m_RandSrc); }
};