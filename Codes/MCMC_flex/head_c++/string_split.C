#include <string>
#include <algorithm>
#include <vector>
#include<iostream>
using namespace std;

void Split(const string& str,
                      vector<string>& tokens,
//                      const string& delimiters = " ")
                      const  string& delimiters )
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
/*
int main()
{
    vector<string> tokens;
    vector<string> v;
    string str("Split me up! Word1 Word2 Word3.");
//unsigned int tokenCount = strtok.countTokens();
    Split(str, tokens,"!");
    v=tokens;
//    copy(tokens.begin(), tokens.end(),v.begin());
    for(int i=0;i<v.size();i++)cout<< v[i]<<" "<<endl;
}

*/
