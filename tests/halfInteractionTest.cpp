/*
 * halfInteractionTest.cpp
 *
 *  Created on: 15 May 2020
 *  Author: Richard Thompson (ithompson@hbku.edu.qa)
 */

#include "hicupData.h"
#include "BinomData.h"
#include "Utils.h"
#include <vector>
#define BOOST_TEST_MODULE SystemCheck
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(SystemCheck)

BOOST_AUTO_TEST_CASE(systemCheck) {

	halfInteraction I = halfInteraction("chr1", 1357908642);
	halfInteraction T(I);

	BOOST_CHECK_EQUAL(I.getInt(), "chr1:1357908642");
	BOOST_CHECK_EQUAL(T.getInt(), "chr1:1357908642");

}

BOOST_AUTO_TEST_SUITE_END()




