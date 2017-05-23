// Mina Jafrai
// 12-14-2015

#include "BindingSiteClass.h"
#include <string>

BindingSiteClass::BindingSiteClass()
{
    for (int i=0; i<3; i++)
    {
        mCoordinate[i] = 0.0;
    }
}

BindingSiteClass::BindingSiteClass(std::string inType, double inX, double inY, double inZ):
    mType("TT")
{
    for (int i=0; i<3; i++)
    {
        mCoordinate[i] = 0.0;
    }
    setType(inType);
    setCoordinates(inX, inY, inZ);
}

std::string BindingSiteClass::getType() const
{
    return (mType);
}

double BindingSiteClass::getX() const
{
    return (mCoordinate[0]);
}

double BindingSiteClass::getY() const
{
    return (mCoordinate[1]);
}

double BindingSiteClass::getZ() const
{
    return (mCoordinate[2]);
}

void BindingSiteClass::getXYZ(double* carts)
{
    carts[0] = mCoordinate[0];
    carts[1] = mCoordinate[1];
    carts[2] = mCoordinate[2];
}

void BindingSiteClass::setType (std::string inType)
{
    mType = inType;
}

void BindingSiteClass::setCoordinates(double inX, double inY, double inZ)
{
    mCoordinate[0] = inX;
    mCoordinate[1] = inY;
    mCoordinate[2] = inZ;
}
