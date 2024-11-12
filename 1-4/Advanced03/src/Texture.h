#pragma once

#include "glm/glm.hpp"
#include <vector>

class Texture
{
public:
	static void ClearTextureCache()
	{
		for (int i=0; i<(int)s_TextureCache.size(); i++)
			delete s_TextureCache[i];
		s_TextureCache.clear();
	}

	enum Texture_Type
	{
		Checker_Texture,
		Solid_Checker_Texture,
		Noise_Texture,
		Marble_Texture,
		Wooden_Texture,
		Image_Texture,
		Turbulence_Texture,
		Dturbulence_Texture,
	};

	Texture()
	{
	}

	Texture(const Texture *t)
	{
	}

	Texture(const Texture &t)
	{
	}

	virtual ~Texture()
	{
	}

	virtual glm::vec3 getColor(const glm::vec3 &uvw) const
	{
		return glm::vec3(0.f, 0.f, 0.f);
	}

	inline glm::vec3 getColor(float u, float v, float w) const
	{
		return getColor(glm::vec3(u,v,w));
	}

	inline glm::vec3 getColor(const glm::vec2 &uv) const
	{
		return getColor(uv.x, uv.y, 0.f);
	}

	virtual Texture_Type getTextureType() const = 0;
	virtual bool isSolidTexture() const = 0;
	
protected:
	static std::vector<Texture *> s_TextureCache;

};
