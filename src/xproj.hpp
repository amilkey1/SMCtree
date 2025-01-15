#pragma once

namespace proj {

    class XProj : public exception {
        public:
                                XProj() throw() {}
                                XProj(const string s) throw() : _msg() {_msg = s;}
                                XProj(const format & f) throw() : _msg() {_msg = str(f);}
            virtual             ~XProj() throw() {}
            const char *        what() const throw() {return _msg.c_str();}

        private:

            string         _msg;
    };

}
