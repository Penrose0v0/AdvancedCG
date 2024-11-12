#pragma once

#include "HitRecord.h"
#include <vector>

class Material
{
protected:
	typedef glm::vec3 vec3;
	typedef float Real;

	Material()
		: ambient(0,0,0)
	{
	}

	Material(const glm::vec3 &_ambient)
		: ambient(_ambient)
	{
	}

	Material(const Material *m)
		: ambient(m->ambient)
	{
	}

	Material(const Material &m)
		: ambient(m.ambient)
	{
	}

public:
	enum Material_Type
	{
		Pseudo_Normal_Color_Type,
		Ambient_Type,
		Diffuse_Type,
		Blinn_Phong_Type,
		Textured_Type,
		Perfect_Specular_Type,
		Specular_Refraction_Type
	};

	static void ClearMaterialCache()
	{
		for (int i=0; i<(int)material_cache.size(); i++)
			delete material_cache[i];
		material_cache.clear();
	}

	static Material *CreateMaterial()
	{
		Material *m = new Material();
		material_cache.push_back( m );
		return m;
	}

	static Material *CreateMaterial(const vec3 &_ambient)
	{
		Material *m = new Material(_ambient);
		material_cache.push_back( m );
		return m;
	}

	static Material *CreateMaterial(const Material *_m)
	{
		Material *m = new Material(_m);
		material_cache.push_back( m );
		return m;
	}

	static Material *CreateMaterial(const Material &_m)
	{
		Material *m = new Material(_m);
		material_cache.push_back( m );
		return m;
	}

	virtual ~Material()
	{
	}

	virtual Material_Type getMaterialType() const { return Ambient_Type; }

	virtual void clone(Material **m) const { *m = new Material(*this); }

	vec3 getAmbientCoeff() const { return ambient; }
	void setAmbientCoeff(const vec3 &a) { ambient = a; }
	void setAmbientCoeff(Real r, Real g, Real b) { ambient = glm::vec3(r,g,b); }

protected:
	vec3 ambient;

	static std::vector<Material *> material_cache;

};

