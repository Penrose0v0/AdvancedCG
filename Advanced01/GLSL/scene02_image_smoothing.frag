#version 130

in vec2 outTexCoord;
out vec4 fragColor;

// TODO: uncomment these lines
uniform sampler2D tex;
uniform int halfKernelSize;
uniform float uScale;
uniform float vScale;

void main()
{
	// TODO: write an appropriate code here
	for (int dx = -halfKernelSize; dx <= halfKernelSize; dx++) {
		for (int dy = -halfKernelSize; dy <= halfKernelSize; dy++) {
			vec2 offset = vec2(dx * uScale, dy * vScale); 
			vec4 currentColor = texture2D(tex, outTexCoord + offset); 
			fragColor += currentColor; 
		}
	}
	fragColor /= fragColor.a; 

}
