#pragma once

#include <cstdio>
#include <vector>
#include <string>

class PathFinder
{
public:
	void addSearchPath(const std::string& dirname) { m_SearchDirs.push_back(dirname); }

	std::string find(const std::string& filename, const std::string& separator = "/")
	{
		FILE* fp = 0;
		std::string path_name;

		for (unsigned int i = 0; i < m_SearchDirs.size(); i++)
		{
			path_name = m_SearchDirs[i] + separator + filename;
#ifdef _WIN32
			::fopen_s(&fp, path_name.c_str(), "r");
#else
			fp = fopen(path_name.c_str(), "r");
#endif

			if (fp != 0)
			{
				fclose(fp);
				return path_name;
			}
		}

		path_name = filename;
#ifdef _WIN32
		::fopen_s(&fp, path_name.c_str(), "r");
#else
		fp = fopen(path_name.c_str(), "r");
#endif

		if (fp != 0)
		{
			fclose(fp);
			return path_name;
		}

		//printf("file not found: %s\n", filename.c_str());

		return "";
	}

private:
	std::vector<std::string> m_SearchDirs;
};