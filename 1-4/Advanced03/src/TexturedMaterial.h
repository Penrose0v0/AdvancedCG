#pragma once

#include "Material.h"
#include "Texture.h"
#include <algorithm>

// Similar to PhongMaterial, but the diffuse component is replaced by texture color

class TexturedMaterial : public Material
{
private:
	TexturedMaterial()
		: shininess(32), ks(1,1,1)
	{
	}

	TexturedMaterial(const TexturedMaterial &m)
		: shininess(m.shininess), ks(m.ks), Material(m)
	{
		texture = m.texture;
	}

public:
	static TexturedMaterial *CreateMaterial()
	{
		TexturedMaterial *m = new TexturedMaterial();
		material_cache.push_back( m );
		return m;
	}

	static TexturedMaterial *CreateMaterial(const TexturedMaterial *_m)
	{
		TexturedMaterial *m = new TexturedMaterial(*_m);
		material_cache.push_back( m );
		return m;
	}

	~TexturedMaterial()
	{
		if (texture) delete texture;
	}

	Material_Type getMaterialType() const { return Textured_Type; }

	void setTexture(Texture *t) { texture = t; }
	const Texture *getTexture() const { return texture; }

	Real getShininess() const { return shininess; }
	void setShininess(Real s) { shininess = s; }

	vec3 getPhongCoeff() const { return ks; }
	void setPhongCoeff(const vec3 &k) { ks = k; }
	void setPhongCoeff(Real r, Real g, Real b) { ks = vec3(r,g,b); }

private:
	Real shininess; // power of cosine lobe
	vec3 ks;	// phong
	Texture *texture;	// pointer for texture, supplies the diffuse color

};
