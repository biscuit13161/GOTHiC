/*
 * unitTest1.cpp
 *
 *  Created on: 11 May 2020
 *      Author: rich
 */

[TestClass]
public class HiCUP_Tests
{
    [TestMethod]
    public void TestMethod1()
    {
        vector<Interaction> interactions;
        interactions.push_back(Interactions("chr1","chr2",12553,15273));
        interactions.push_back(Interactions("chr1","chr2",12553,15273));
        interactions.push_back(Interactions("chr1","chr2",12553,15273));
        interactions.push_back(Interactions("chr6","chrX",125585523,1063441));
        interactions.push_back(Interactions("chr6","chrX",125585523,1063441));
        interactions.push_back(Interactions("chr10","chr5",1064473,1505273));

        vector<Interaction> results;
        results.push_back(Interactions("chr1","chr2",12553,15273,3));
        results.push_back(Interactions("chr6","chrX",125585523,1063441,2));
        results.push_back(Interactions("chr10","chr5",1064473,1505273,1));

        countDuplicates(vector<Interaction> & interactions);

        Assert.AreEqual(results, interactions);
    }
}


