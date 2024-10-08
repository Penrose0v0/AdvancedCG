#version 130

in vec2 outTexCoord;
out vec4 fragColor;

// TODO: uncomment these lines
uniform vec4 checkerColor0;
uniform vec4 checkerColor1;
uniform vec2 checkerScale;

void main()
{
	// TODO: write an appropriate code here
	vec2 checkerIndex = floor(outTexCoord / (checkerScale / 2));  // Get index
	float checker = mod(checkerIndex.x + checkerIndex.y, 2.0);  // Choose color
    fragColor = (checker < 1.0) ? checkerColor0 : checkerColor1;  // Set color
}
