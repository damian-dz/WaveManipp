#pragma once

#ifndef WAVE_BUILDER_H
#define WAVE_BUILDER_H

#include "Wave.hpp"

namespace wm {
/*!
 * \brief A class that is useful for joining many Wave objects together.
 * \details This class makes it possible to join multiple Wave objects without memory reallocation overhead.
 */
class WaveBuilder
{
public:
    WaveBuilder();
};
}

#endif // !WAVE_BUILDER_H
