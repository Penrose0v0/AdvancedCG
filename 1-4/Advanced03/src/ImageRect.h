#pragma once

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>

template <class T>
class ImageRect
{
protected:
	int width;
	int height;
	T *data;

public:
	ImageRect()
		: data(0)
	{
	}

	ImageRect(int w, int h)
		: width(w), height(h)
	{
		data = new T[ width*height ];
	}

	ImageRect(int w, int h, const T &v)
		: width(w), height(h)
	{
		data = new T[ width*height ];
		fill( v );
	}

	ImageRect(int w, int h, T *_data)
		: width(w), height(h)
	{
		assert(_data);

		data = new T[ width*height ];
		memcpy(data, _data, width*height*sizeof(T));
	}

	ImageRect(const ImageRect &img)
		: width(img.width), height(img.height)
	{
		data = new T[ width*height ];

		if (img.data)
		{
			memcpy(data, img.data, width*height*sizeof(T));
		}
	}

	~ImageRect()
	{
		if (data) delete [] data;
	}

	inline int getWidth() const { return width; }
	inline int getHeight() const { return height; }
	inline int getImageSize() const { return width*height; }
	inline T * getData() { return data; }
	inline T const * getData() const { return data; }

	inline T * getScanline(int yi) { return data + yi*width; }
	inline T const * getScanline(int yi) const { return data + yi*width; }

	void allocate(int w, int h)
	{
		if (data) delete [] data;
		data = new T[w*h];
		width = w;
		height = h;
	}

	void allocate(int w, int h, T *_data)
	{
		if (data) delete [] data;
		data = new T[w*h];
		memcpy(data, _data, w*h*sizeof(T));
		width = w;
		height = h;
	}

	void copy(const ImageRect &img)
	{
		if (data) delete [] data;

		width = img.width;
		height = img.height;
		data = new T[width*height];
		memcpy(data, img.data, width*height*sizeof(T));
	}

	inline const T getValue(int xi, int yi) const
	{
		assert(data);
		assert(0 <= xi && xi < width);
		assert(0 <= yi && yi < height);

		return data[xi + width*yi];
	}

	const T bilinearInterp(float x, float y) const
	{
		assert(data);

		// clamped to edge

		const float _x = std::max(0.f, std::min(x, (float)(width-1)));
		const float _y = std::max(0.f, std::min(y, (float)(height-1)));

		const int x0 = (int)floorf(_x);
		const int x1 = std::min(x0+1, width-1);
		const float s = _x - x0;

		const int y0 = (int)floorf(_y);
		const int y1 = std::min(y0+1, height-1);
		const float t = _y - y0;

		const T v00 = data[x0 + width*y0];
		const T v10 = data[x1 + width*y0];
		const T v01 = data[x0 + width*y1];
		const T v11 = data[x1 + width*y1];

		const T r0 = v00 + s*(v10 - v00);
		const T r1 = v01 + s*(v11 - v01);

		return r0 + t*(r1 - r0);
	}

	inline void setValue(int xi, int yi, T& v)
	{
		assert(data);
		assert(0 <= xi && xi < width);
		assert(0 <= yi && yi < height);

		data[xi + width*yi] = v;
	}

	inline T& operator()(int xi, int yi)
	{
		assert(data);
		assert(0 <= xi && xi < width);
		assert(0 <= yi && yi < height);

		return data[xi + width*yi];
	}

	inline const T& operator()(int xi, int yi) const
	{
		assert(data);
		assert(0 <= xi && xi < width);
		assert(0 <= yi && yi < height);

		return data[xi + width*yi];
	}

	inline ImageRect& operator=(const ImageRect &img)
	{
		delete [] data;
		width = img.width;
		height = img.height;
		data = new T[width*height];
		memcpy(data, img.data, width*height*sizeof(T));
		return *this;
	}

	inline ImageRect& operator+=(const float f)
	{
		const int size = width*height;
		for (int i=0; i<size; i++)
			data[i] += f;
		return *this;
	}

	inline ImageRect& operator-=(const float f)
	{
		const int size = width*height;
		for (int i=0; i<size; i++)
			data[i] -= f;
		return *this;
	}

	inline ImageRect& operator*=(const float f)
	{
		const int size = width*height;
		for (int i=0; i<size; i++)
			data[i] *= f;
		return *this;
	}

	inline ImageRect& operator/=(const float f)
	{
		const int size = width*height;
		const float fInv = 1.f/f;
		for (int i=0; i<size; i++)
			data[i] *= fInv;
		return *this;
	}

	inline ImageRect& operator+=(const ImageRect &img)
	{
		assert(width == img.width);
		assert(height == img.height);
		const int size = width*height;
		for (int i=0; i<size; i++)
			data[i] += img.data[i];
		return *this;
	}

	inline ImageRect& operator-=(const ImageRect &img)
	{
		assert(width == img.width);
		assert(height == img.height);
		const int size = width*height;
		for (int i=0; i<size; i++)
			data[i] -= img.data[i];
		return *this;
	}

	inline ImageRect& operator*=(const ImageRect &img)
	{
		assert(width == img.width);
		assert(height == img.height);
		const int size = width*height;
		for (int i=0; i<size; i++)
			data[i] *= img.data[i];
		return *this;
	}

	inline ImageRect& operator/=(const ImageRect &img)
	{
		assert(width == img.width);
		assert(height == img.height);
		const int size = width*height;
		for (int i=0; i<size; i++)
			data[i] /= img.data[i];
		return *this;
	}

	void getSubImage(int x0, int y0, int w, int h, ImageRect &img)
	{
		assert(x0+w < width);
		assert(y0+h < height);
		assert(width);
		assert(height);
		assert(w);
		assert(h);

		if (img.data) delete [] img.data;
		img.data = new T[w*h];
		img.width = w;
		img.height = h;

		for (int yi=0; yi<h; yi++)
		{
			memcpy(&img.data[w*yi], &data[x0+width*(y0+yi)], w*sizeof(T));
		}
	}

	void setSubImage(int x0, int y0, ImageRect &img)
	{
		assert(data);
		assert(x0+img.width < width);
		assert(y0+img.height < height);
		assert(width);
		assert(height);
		assert(img.width);
		assert(img.height);

		for (int yi=0; yi<img.height; yi++)
		{
			memcpy(&data[x0+width*(y0+yi)], &img.data[img.width*yi], img.width*sizeof(T));
		}
	}

	void getRegionData(int x0, int y0, int w, int h, T * _data)
	{
		assert(data);
		assert(x0+w < width);
		assert(y0+h < height);
		assert(width);
		assert(height);
		assert(w);
		assert(h);

		if (_data) delete [] _data;
		_data = new T[w*h];

		for (int yi=0; yi<h; yi++)
		{
			memcpy(&_data[w*yi], &data[x0+width*(y0+yi)], w*sizeof(T));
		}
	}

	void setRegionData(int x0, int y0, int w, int h, T * _data)
	{
		assert(data);
		assert(x0+w < width);
		assert(y0+h < height);
		assert(width);
		assert(height);
		assert(w);
		assert(h);

		for (int yi=0; yi<h; yi++)
		{
			memcpy(&data[x0+width*(y0+yi)], &_data[w*yi], w*sizeof(T));
		}
	}

	void fill(const T& v)
	{
		assert(data);

		for (int yi=0; yi<height; yi++)
			for (int xi=0; xi<width; xi++)
				data[xi + width*yi] = v;
	}

	void trim(int x0, int y0, int w, int h)
	{
		assert(data);
		assert(x0+w < width);
		assert(y0+h < height);
		assert(width);
		assert(height);
		assert(w);
		assert(h);

		T * _tmp = new T[w*h];

		for (int yi=0; yi<h; yi++)
		{
			memcpy(&_tmp[w*yi], &data[x0 + width*(y0+yi)], w*sizeof(T));
		}

		delete [] data;
		data = new T[w*h];
		memcpy(data, _tmp, w*h*sizeof(T));
		width = w;
		height = h;

		delete [] _tmp;
	}

	
};


