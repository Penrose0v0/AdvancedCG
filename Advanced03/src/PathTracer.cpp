#include <GL/glew.h>
#include "PathTracer.h"
#include "GLFW/glfw3.h"
#include "arcball_camera.h"
#include "Scene.h"
#include "HitRecord.h"
#include <iostream>
#include <chrono>

#include "DiffuseMaterial.h"
#include "BlinnPhongMaterial.h"
#include "PerfectSpecularMaterial.h"
#include "SpecularRefractionMaterial.h"

using namespace std;
using namespace glm;

float PathTracer::s_Gamma = 2.2f;
int PathTracer::s_MaxRecursionDepth = 32;
int PathTracer::s_MinRecursionDepth = 5;
int PathTracer::s_NumSamplesPerPixel = 100;
int PathTracer::s_NumSamplesPerUpdate = 5;

extern GLFWwindow* g_pWindow;
extern ArcballCamera g_Camera;
extern Scene g_Scene;
extern int g_WindowWidth, g_WindowHeight;
extern mat4 g_ProjMatrix;
extern bool g_KeepTracing;

void PathTracer::renderScene()
{
	const auto tStart = chrono::system_clock::now();

	vec3 xAxis, yAxis, zAxis, eye;
	g_Camera.getEyeCoordinateSystem(xAxis, yAxis, zAxis, eye);

	m_FrameBuffer.allocate(g_WindowWidth, g_WindowHeight);

	const float halfWidth = 0.5f * g_WindowWidth;
	const float halfHeight = 0.5f * g_WindowHeight;
	const float screenDist = halfHeight * g_ProjMatrix[1][1];

	int nRemainingSamples = s_NumSamplesPerPixel;

	while (nRemainingSamples > 0)
	{
		const int nSamplesDone = s_NumSamplesPerPixel - nRemainingSamples;
		const int nNewRemainingSamples = std::max(nRemainingSamples - s_NumSamplesPerUpdate, 0);
		const int nNewSamples = nRemainingSamples - nNewRemainingSamples;
		const int nNewSamplesDone = nSamplesDone + nNewSamples;

#pragma omp parallel for
		for (int yi = 0; yi < g_WindowHeight; ++yi)
		{
			for (int xi = 0; xi < g_WindowWidth; ++xi)
			{
				vec3 pixelColor = float(nSamplesDone) * m_FrameBuffer(xi, yi);

				for (int si = 0; si < nNewSamples; ++si)
				{
					const vec3 dir = (xi + frand() - halfWidth) * xAxis + (yi + frand() - halfHeight) * yAxis - screenDist * zAxis;
					pixelColor += traceRec(Ray(eye, glm::normalize(dir)), 0);
				}

				m_FrameBuffer(xi, yi) = pixelColor / float(nNewSamplesDone);
			}
		}

		nRemainingSamples = nNewRemainingSamples;

		renderIntermediateFrame();

		const auto tNow = chrono::system_clock::now();
		const auto elapsed = chrono::duration_cast<chrono::milliseconds>(tNow - tStart).count() / 1000.f;
		cerr << __FUNCTION__ << ": " << nNewSamplesDone << "/" << s_NumSamplesPerPixel << " samples (" << elapsed << " sec)" << endl;
	}
}

void PathTracer::renderFrame()
{
	if (!m_FrameBufferTexID) return;

	if (!m_pGammaShader)
		initShader();

	const float uOffset = (m_isNVIDIADriver) ? -0.5f / g_WindowWidth : 0.f;
	const float vOffset = (m_isNVIDIADriver) ? -0.5f / g_WindowHeight : 0.f;

	glBindTexture(GL_TEXTURE_2D, m_FrameBufferTexID);

	m_pGammaShader->use();
	m_pGammaShader->sendUniform1i("tex", 0);
	m_pGammaShader->sendUniform1f("gamma", s_Gamma);
	m_pGammaShader->sendUniform2f("offset", uOffset, vOffset);

	glBegin(GL_QUADS);
	glTexCoord2f(0, 0);	glVertex2f(-1, -1);
	glTexCoord2f(1, 0);	glVertex2f( 1, -1);
	glTexCoord2f(1, 1);	glVertex2f( 1,  1);
	glTexCoord2f(0, 1);	glVertex2f(-1,  1);
	glEnd();

	m_pGammaShader->disable();

	glBindTexture(GL_TEXTURE_2D, 0);
}

void PathTracer::initShader()
{
	if (!m_pGammaShader)
		m_pGammaShader = new GLSLProgramObject();

	m_pGammaShader->attachShaderCodeString(
		R"(#version 120
			void main(void)
			{
				gl_TexCoord[0] = gl_MultiTexCoord0;
				gl_Position = gl_Vertex;
			}
		)", GL_VERTEX_SHADER);

	m_pGammaShader->attachShaderCodeString(
		R"(#version 120
			uniform sampler2D tex;
			uniform float gamma;
			uniform vec2 offset;
			void main(void)
			{
				gl_FragColor = vec4(pow(texture2D(tex, gl_TexCoord[0].st + offset).rgb, vec3(1.0/gamma)), 1.0);
			}
		)", GL_FRAGMENT_SHADER);

	m_pGammaShader->link();

	if (!m_pGammaShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": gamma correction shader" << endl;
		m_pGammaShader->printProgramLog();
	}

	// NVIDIA driver requires -0.5 offsets in texture fetch in order to match the path-traced result
	m_isNVIDIADriver = strncmp((const char*)glGetString(GL_VENDOR), "NVIDIA", sizeof("NVIDIA") - 1) == 0;
}

glm::vec3 PathTracer::traceRec(const Ray& ray, int recursionDepth)
{
	if (recursionDepth > s_MaxRecursionDepth)
		return g_Scene.getBackgroundColor(ray);

	const float tEpsilon = 0.01f;
	const float tInfinity = 1.0e+10f;

	HitRecord record;
	record.m_ParamT = tInfinity;

	bool hitObject = false;

	for (int oi = 0; oi < g_Scene.getNumObjects(); oi++)
	{
		GeometricObject* o = g_Scene.getObject(oi);

		HitRecord tmpRec;
		const bool isHit = o->hit(ray, tEpsilon, tInfinity, tmpRec);

		if (isHit && record.m_ParamT > tmpRec.m_ParamT)
		{
			record = tmpRec;
			hitObject = true;
		}
	}

	if (!hitObject)
		return g_Scene.getBackgroundColor(ray);

	const Material::Material_Type matType = record.m_pMaterial->getMaterialType();

	if (matType == Material::Pseudo_Normal_Color_Type)
	{
		return 0.5f * record.m_Normal + vec3(0.5f);
	}
	else if (matType == Material::Diffuse_Type)
	{
		// TODO: implement Lambert (diffuse) reflection with BRDF importance sampling

		// HINT: local coordinate system can be defined using the following function:
		//vec3 xLocal, yLocal, zLocal;
		//calcLocalCoordinateSystem(record.m_Normal, ray.getUnitDir(), xLocal, yLocal, zLocal);
		const float M_PI = 3.14159265358979323;
		vec3 normal = normalize(record.m_Normal);
		vec3 incomingDir = normalize(ray.getUnitDir());

		// Create an orthonormal basis (xLocal, yLocal, zLocal)
		vec3 zLocal = normal;
		vec3 xLocal = normalize(glm::cross((fabs(zLocal.x) > 0.1f ? vec3(0.0f, 1.0f, 0.0f) : vec3(1.0f, 0.0f, 0.0f)), zLocal));
		vec3 yLocal = glm::cross(zLocal, xLocal);

		// Cosine-weighted hemisphere sampling
		float r1 = frand();
		float r2 = frand();
		float phi = 2.0f * M_PI * r1;
		float cosTheta = sqrt(1.0f - r2);
		float sinTheta = sqrt(r2);

		// Convert spherical coordinates to Cartesian coordinates in the local system
		vec3 sampledDir = normalize(
			cos(phi) * sinTheta * xLocal +
			sin(phi) * sinTheta * yLocal +
			cosTheta * zLocal
		);

		// Ensure the sampled direction is in the same hemisphere as the normal
		if (dot(sampledDir, normal) < 0.0f)
			sampledDir = -sampledDir;

		// Offset the hit position slightly to prevent self-intersection
		vec3 newOrigin = record.m_HitPos + normal * 0.001f;

		// Create the new ray
		Ray newRay(newOrigin, sampledDir);

		// Retrieve the diffuse reflectance
		DiffuseMaterial* diffuseMat = static_cast<DiffuseMaterial*>(record.m_pMaterial);
		vec3 diffuseCoeff = diffuseMat->getDiffuseCoeff();

		// Calculate the BRDF value for Lambertian reflection
		// BRDF = diffuseCoeff / PI
		vec3 brdf = diffuseCoeff / M_PI;

		// Calculate the cosine of the angle between the normal and the sampled direction
		float cosAlpha = std::max(dot(sampledDir, normal), 0.0f);

		// Russian Roulette termination
		float probability = 1.0f;
		if (recursionDepth >= s_MinRecursionDepth)
		{
			// Probability based on the maximum component of diffuse reflectance
			probability = std::max(diffuseCoeff.r, std::max(diffuseCoeff.g, diffuseCoeff.b));
			if (frand() >= probability)
				return g_Scene.getBackgroundColor(ray);
		}

		// Trace the new ray recursively
		vec3 incomingRadiance = traceRec(newRay, recursionDepth + 1);

		// Compute the contribution with the BRDF, cosine term, and probability
		vec3 contribution = (brdf * incomingRadiance * cosAlpha) / probability;

		return contribution;
	}
	else if (matType == Material::Blinn_Phong_Type)
	{
		// TODO: implement Blinn-Phong reflection with BRDF importance sampling
		const float M_PI = 3.14159265358979323;
		BlinnPhongMaterial* blinnMat = static_cast<BlinnPhongMaterial*>(record.m_pMaterial);
		vec3 diffuseCoeff = blinnMat->getDiffuseCoeff();
		vec3 specularCoeff = blinnMat->getSpecularCoeff();
		float shininess = blinnMat->getShininess();

		// Total reflectance
		vec3 totalReflectance = diffuseCoeff + specularCoeff;
		float maxReflectance = std::max(totalReflectance.r, std::max(totalReflectance.g, totalReflectance.b));

		// Russian Roulette termination
		float probability = 1.0f;
		if (recursionDepth >= s_MinRecursionDepth)
		{
			probability = maxReflectance;
			if (frand() >= probability)
				return g_Scene.getBackgroundColor(ray);
		}

		// Determine whether to sample diffuse or specular based on their weights
		float diffuseWeight = length(diffuseCoeff);
		float specularWeight = length(specularCoeff);
		float weightSum = diffuseWeight + specularWeight;

		// Avoid division by zero
		float diffuseProbability = (weightSum > 0.0f) ? (diffuseWeight / weightSum) : 0.0f;
		float specularProbability = (weightSum > 0.0f) ? (specularWeight / weightSum) : 0.0f;

		// Randomly choose to sample diffuse or specular
		float sampleChoice = frand();

		vec3 sampledDir;
		vec3 brdf;
		float pdf = 0.0f;

		// Normalize the normal and incoming direction
		vec3 normal = normalize(record.m_Normal);
		vec3 incomingDir = normalize(ray.getUnitDir());

		// Create an orthonormal basis (xLocal, yLocal, zLocal)
		vec3 zLocal = normal;
		vec3 xLocal = normalize(glm::cross((fabs(zLocal.x) > 0.1f ? vec3(0.0f, 1.0f, 0.0f) : vec3(1.0f, 0.0f, 0.0f)), zLocal));
		vec3 yLocal = glm::cross(zLocal, xLocal);

		if (sampleChoice < diffuseProbability)
		{
			// Diffuse sampling (cosine-weighted hemisphere)

			float r1 = frand();
			float r2 = frand();
			float phi = 2.0f * M_PI * r1;
			float cosTheta = sqrt(1.0f - r2);
			float sinTheta = sqrt(r2);

			// Convert spherical coordinates to Cartesian coordinates in the local system
			sampledDir = normalize(
				cos(phi) * sinTheta * xLocal +
				sin(phi) * sinTheta * yLocal +
				cosTheta * zLocal
			);

			// Ensure the sampled direction is in the same hemisphere as the normal
			if (dot(sampledDir, normal) < 0.0f)
				sampledDir = -sampledDir;

			// BRDF for diffuse
			brdf = diffuseCoeff / M_PI;

			// PDF for cosine-weighted sampling
			pdf = dot(sampledDir, normal) / M_PI;
		}
		else
		{
			// Specular sampling (Blinn-Phong)

			// Sample the half-vector based on the Blinn-Phong distribution
			float r1 = frand();
			float r2 = frand();
			float theta = acos(pow(r1, 1.0f / (shininess + 1.0f)));
			float phi = 2.0f * M_PI * r2;

			float sinTheta = sin(theta);
			float cosTheta = cos(theta);

			// Half-vector in local space
			vec3 halfVector = normalize(
				cos(phi) * sinTheta * xLocal +
				sin(phi) * sinTheta * yLocal +
				cosTheta * zLocal
			);

			// Compute the sampled direction by reflecting the incoming direction around the half-vector
			sampledDir = normalize(reflect(-incomingDir, halfVector));

			// Ensure the sampled direction is in the same hemisphere as the normal
			if (dot(sampledDir, normal) < 0.0f)
				sampledDir = -sampledDir;

			// BRDF for Blinn-Phong specular
			float nDotH = std::max(dot(normal, halfVector), 0.0f);
			brdf = specularCoeff * ((shininess + 2.0f) / (2.0f * M_PI)) * pow(nDotH, shininess);

			// PDF for Blinn-Phong sampling
			pdf = ((shininess + 1.0f) * pow(nDotH, shininess)) / (2.0f * M_PI);
		}

		// Avoid division by zero
		if (pdf <= 0.0f)
			return vec3(0.0f);

		// Offset the hit position slightly to prevent self-intersection
		vec3 newOrigin = record.m_HitPos + normal * 0.001f;

		// Create the new ray
		Ray newRay(newOrigin, sampledDir);

		// Calculate the cosine of the angle between the normal and the sampled direction
		float cosAlpha = std:: max(dot(sampledDir, normal), 0.0f);

		// Trace the new ray recursively
		vec3 incomingRadiance = traceRec(newRay, recursionDepth + 1);

		// Compute the contribution with the BRDF, cosine term, PDF, and probability
		vec3 contribution = (brdf * incomingRadiance * cosAlpha) / (pdf * probability);

		return contribution;

		// return vec3(0.f);
	}
	else if (matType == Material::Perfect_Specular_Type)
	{
		const vec3& specularCoeff = ((PerfectSpecularMaterial*)record.m_pMaterial)->getSpecularCoeff();
		const float russianRouletteProbability = (recursionDepth > s_MinRecursionDepth) ? std::max(specularCoeff.x, std::max(specularCoeff.y, specularCoeff.z)) : 1.f;

		if (frand() >= russianRouletteProbability)
			return g_Scene.getBackgroundColor(ray);

		const vec3 reflectDir = normalize(reflect(ray.getUnitDir(), record.m_Normal));
	
		const vec3 incomingRadiance = traceRec(Ray(record.m_HitPos, reflectDir), recursionDepth + 1);
		const vec3 weight = specularCoeff / russianRouletteProbability;

		return weight * incomingRadiance;
	}
	else if (matType == Material::Specular_Refraction_Type)
	{
		const vec3& specularCoeff = ((SpecularRefractionMaterial*)record.m_pMaterial)->getSpecularCoeff();
		const float russianRouletteProbability = (recursionDepth > s_MinRecursionDepth) ? std::max(specularCoeff.x, std::max(specularCoeff.y, specularCoeff.z)) : 1.f;

		if (frand() >= russianRouletteProbability)
			return g_Scene.getBackgroundColor(ray);

		const float _dot = dot(ray.getUnitDir(), record.m_Normal);

		// is the ray entering or outgoing the ball?
		const bool isEntering = _dot < 0.f;

		const float eta = ((SpecularRefractionMaterial*)record.m_pMaterial)->getRefractionIndex();
		const float relativeIndex = isEntering ? 1 / eta : eta;

		// Schlick's Fresnel approximation

		const float R0 = ((eta - 1.f) * (eta - 1.f)) / ((eta + 1.f) * (eta + 1.f));

		const float c = 1.f - fabsf(_dot);
		const float c2 = c * c;
		const float Re = R0 + (1 - R0) * c2 * c2 * c;
		const float Tr = (1.f - Re) * relativeIndex * relativeIndex;

		const vec3 normal = isEntering ? record.m_Normal : -record.m_Normal;

		const vec3 refractVec = refract(ray.getUnitDir(), normal, relativeIndex);
		const vec3 reflectVec = reflect(ray.getUnitDir(), normal);

		if (refractVec == vec3(0.f)) // total reflection
		{
			const vec3 incomingRadiance = traceRec(Ray(record.m_HitPos, reflectVec), recursionDepth + 1);
			const vec3 weight = specularCoeff / russianRouletteProbability;

			return weight * incomingRadiance;
		}

		if (recursionDepth <= 2)
		{
			const vec3 incomingRadiance = Re * traceRec(Ray(record.m_HitPos, reflectVec), recursionDepth + 1)
										+ Tr * traceRec(Ray(record.m_HitPos, refractVec), recursionDepth + 1);

			const vec3 weight = specularCoeff / russianRouletteProbability;

			return weight * incomingRadiance;
		}
		else
		{
			// apply russian roulette for Fresnel reflection

			const float reflectionProbability = Re;

			if (frand() < reflectionProbability)
			{
				const vec3 incomingRadiance = Re * traceRec(Ray(record.m_HitPos, reflectVec), recursionDepth + 1);
				const vec3 weight = specularCoeff / (reflectionProbability * russianRouletteProbability);

				return weight * incomingRadiance;
			}
			else
			{
				const vec3 incomingRadiance = Tr * traceRec(Ray(record.m_HitPos, refractVec), recursionDepth + 1);
				const vec3 weight = specularCoeff / ((1.f - reflectionProbability) * russianRouletteProbability);

				return weight * incomingRadiance;
			}
		}
	}

	return vec3(0.f);
}

void PathTracer::updateFrameBufferTexture()
{
	if (!m_FrameBufferTexID) glGenTextures(1, &m_FrameBufferTexID);
	glBindTexture(GL_TEXTURE_2D, m_FrameBufferTexID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, g_WindowWidth, g_WindowHeight, 0, GL_RGB, GL_FLOAT, m_FrameBuffer.getData());
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void PathTracer::renderIntermediateFrame()
{
	updateFrameBufferTexture();

	if (!g_KeepTracing)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		renderFrame();
		glfwSwapBuffers(g_pWindow);
	}
}

void PathTracer::calcLocalCoordinateSystem(const vec3& normal, const vec3& inDir, vec3& xLocal, vec3& yLocal, vec3& zLocal) const
{
#if 0
	yLocal = dot(normal, inDir) < 0.f ? normal : -normal;
	xLocal = normalize(cross(fabsf(normal.x) > 0.001f ? vec3(0, 1, 0) : vec3(1, 0, 0), normal));
	zLocal = cross(xLocal, yLocal);
#else
	yLocal = normal;

	if (fabsf(glm::dot(normal, inDir)) < 0.9f)
		xLocal = glm::normalize(glm::cross(inDir, normal));
	else
		xLocal = glm::normalize((fabsf(normal.y) > 0.1f) ? glm::vec3(normal.y, -normal.x, 0.f) : glm::vec3(normal.z, 0.f, -normal.x));

	zLocal = glm::cross(xLocal, yLocal);
#endif
}
