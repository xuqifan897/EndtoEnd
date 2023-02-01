#ifndef __CONFIGPARSER_H__
#define __CONFIGPARSER_H__

#include <string>
#include "rapidjson/document.h"
#include <iostream>
#include <vector>

#include "DoseCalcIO/io_helpers.h"

// provides a collection of key-value strings extracted from a json formatted config file
rapidjson::Document read_config(const std::string& fname, int verbose=false);

template <typename T>
int static get_and_assign_scalar(T& assign, rapidjson::Document& doc, const std::string& key, int verbose=false) {
    int status = 0;
    rapidjson::Value::ConstMemberIterator itr = doc.FindMember(key.c_str());
    if (itr != doc.MemberEnd()) {
        if (itr->value.IsNumber()) {
            assign = itr->value.Get<T>();
            doc.RemoveMember(key.c_str());
            status = true;
        }
    }

    if (verbose) {
        if (status) { std::cout << "- FOUND   \""<<key<<"\" = " << assign << std::endl; }
        else { std::cout << "- MISSING \""<<key<<"\"" << std::endl; }
    }
    return status;
}
template <> int get_and_assign_scalar<std::string>(std::string& assign, rapidjson::Document& doc, const std::string& key, int verbose) {
    int status = 0;
    rapidjson::Value::ConstMemberIterator itr = doc.FindMember(key.c_str());
    if (itr != doc.MemberEnd()) {
        if (itr->value.IsString()) {
            assign = std::string(itr->value.GetString());
            doc.RemoveMember(key.c_str());
            status = true;
        }
    }

    if (verbose) {
        if (status) { std::cout << "- FOUND   \""<<key<<"\" = " << assign << std::endl; }
        else { std::cout << "- MISSING \""<<key<<"\"" << std::endl; }
    }
    return status;
}
template <> int get_and_assign_scalar<int>(int& assign, rapidjson::Document& doc, const std::string& key, int verbose) {
    int status = 0;
    rapidjson::Value::ConstMemberIterator itr = doc.FindMember(key.c_str());
    if (itr != doc.MemberEnd()) {
        if (itr->value.IsString()) {
            assign = dcio::tolower(std::string(itr->value.GetString())).compare("true") == 0 ?  1 : 0;
        } else if (itr->value.IsBool()) {
            assign = (int)itr->value.GetBool();
        } else if (itr->value.IsNumber()) {
            assign = (int)itr->value.GetFloat();
        }
        doc.RemoveMember(key.c_str());
        status = true;
    }

    if (verbose) {
        if (status) { std::cout << "- FOUND   \""<<key<<"\" = " << assign << std::endl; }
        else { std::cout << "- MISSING \""<<key<<"\"" << std::endl; }
    }
    return status;
}
template <> int get_and_assign_scalar<bool>(bool& assign, rapidjson::Document& doc, const std::string& key, int verbose) {
    int status = 0;
    rapidjson::Value::ConstMemberIterator itr = doc.FindMember(key.c_str());
    if (itr != doc.MemberEnd()) {
        if (itr->value.IsString()) {
            assign = dcio::tolower(std::string(itr->value.GetString())).compare("true") == 0 ?  true : false;
        } else if (itr->value.IsBool()) {
            assign = itr->value.GetBool();
        } else if (itr->value.IsNumber()) {
            assign = itr->value.GetFloat() == 0.f ? true : false;
        }
        doc.RemoveMember(key.c_str());
        status = true;
    }

    if (verbose) {
        if (status) { std::cout << "- FOUND   \""<<key<<"\" = " << assign << std::endl; }
        else { std::cout << "- MISSING \""<<key<<"\"" << std::endl; }
    }
    return status;
}

template<typename baseT, typename vectT, unsigned int N>
static int get_and_assign_vector(vectT& assign, rapidjson::Document& doc, const std::string& key, int verbose=false) {
    int status = 0;
    rapidjson::Value::ConstMemberIterator itr = doc.FindMember(key.c_str());
    if ( itr != doc.MemberEnd() ) {
        if (itr->value.IsNumber()) {
            baseT val = itr->value.GetDouble();
            assign = vectT{}; // empty existing value
            for (int i=0; i<N; ++i) {
                *(((baseT*)&assign)+i) = static_cast<baseT>(val);
            }
            status = true;
        }
        else if (itr->value.IsArray() && itr->value.Size()>=2) {
            auto arr = itr->value.GetArray();
            assign = vectT{};
            for (int i=0; i<N; ++i) {
                *(((baseT*)&assign)+i) = static_cast<baseT>(arr[i].GetDouble());
            }
            status = true;
        } else { std::cout << "invalid type/format for member \""<<key<<"\"\n"; }
    }

    if (verbose) {
        if (status) {
            std::cout << "- FOUND   \""<<key<<"\" = {";
            for (int i=0; i<N; ++i) {
                if (i != 0) { std::cout << ", "; }
                std::cout << *(((baseT*)&assign)+i);
            }
            std::cout << "}" << std::endl;
        } else {
            std::cout << "- MISSING \""<<key<<"\"" << std::endl;
        }
    }
    return status;
}

template<unsigned int N>
static int get_and_assign_vector_string(std::vector<std::string>& assign, rapidjson::Document& doc, const std::string& key, int verbose) {
    int status = 0;
    rapidjson::Value::ConstMemberIterator itr = doc.FindMember(key.c_str());
    if ( itr != doc.MemberEnd() ) {
        assign = std::vector<std::string>(N); // empty existing value
        if (itr->value.IsString()) {
            std::string val = itr->value.GetString();
            for (int i=0; i<N; ++i) {
                assign[i] = std::string(val);
            }
            status = true;
        }
        else if (itr->value.IsArray() && itr->value.Size()>=2) {
            auto arr = itr->value.GetArray();
            for (int i=0; i<N; ++i) {
                assign[i] = std::string(arr[i].GetString());
            }
            status = true;
        } else { std::cout << "invalid type/format for member \""<<key<<"\"\n"; }
    }

    if (verbose) {
        if (status) {
            std::cout << "- FOUND   \""<<key<<"\" = {";
            for (int i=0; i<N; ++i) {
                if (i != 0) { std::cout << ", "; }
                std::cout << assign[i];
            }
            std::cout << "}" << std::endl;
        } else {
            std::cout << "- MISSING \""<<key<<"\"" << std::endl;
        }
    }
    return status;
}

#endif // __CONFIGPARSER_H__
