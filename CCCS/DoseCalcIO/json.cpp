#include "json.h"

#include <cstdio>
#include <exception>

// json parser (namespace rapidjson::)
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

typedef std::vector<std::string> strvect;

strvect getROINamesFromJSON(std::string json_path) {
    // parse json
    strvect roi_names;

    FILE* fp;
    if ( (fp = fopen(json_path.c_str(), "r")) != NULL ) {
        char readbuf[5000];
        rapidjson::FileReadStream is{fp, readbuf, sizeof(readbuf)};
        rapidjson::Document jdoc;
        jdoc.ParseStream(is);

        if (jdoc.IsObject()) {
            rapidjson::Value::ConstMemberIterator ptv_itr = jdoc.FindMember("ptv");
            if (ptv_itr != jdoc.MemberEnd() && ptv_itr->value.IsString()) {
                roi_names.emplace_back(ptv_itr->value.GetString());
            } else {
                throw std::exception();
            }

            rapidjson::Value::ConstMemberIterator oar_itr = jdoc.FindMember("oar");
            if (oar_itr != jdoc.MemberEnd()) {
                if (oar_itr->value.IsArray()) {
                    for (auto& v : oar_itr->value.GetArray()) {
                        if (v.IsString()) {
                            roi_names.emplace_back( v.GetString() );
                        }
                    }
                } else if (oar_itr->value.IsString()) {
                    roi_names.emplace_back( oar_itr->value.GetString() );
                }
            }
        }
        fclose(fp);
    }
    return roi_names;
}

