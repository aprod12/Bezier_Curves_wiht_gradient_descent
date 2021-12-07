#pragma once
enum Dimension
{
	TWO_DIMENSION = 2,
	THREE_DIMENSION = 3
};


class IDrawable
{
public:
	virtual size_t getDegree() = 0;
};

