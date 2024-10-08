#version 130

// TODO: define "in" variables

out vec4 fragColor;

uniform sampler2DShadow shadowTex;
uniform vec2 texMapScale;

// TODO: define uniform variables

float offsetLookup(sampler2DShadow map, vec4 loc, vec2 offset)
{
	return textureProj(map, vec4(loc.xy + offset * texMapScale * loc.w, loc.z, loc.w));
}

float samplingLimit = 3.5;

void main()
{
	// HINT: The visibility (i.e., shadowed or not) can be fetched using the "offsetLookup" function
	//       in the double loops. The y and x loops ranges from -samplingLimit to samplingLimit
	//       with a stepsize 1.0. Finally, the summed visibility is divided by 64 to normalize.

	// TODO: rewirte this function
	fragColor = vec4(0, 1, 0, 1);
}
