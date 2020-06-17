/*
 *  AUTHOR
 *	Richard Thompson, ithompson@hbku.edu.qa
 *	October 23, 2020.
 *
 *	Copyright (C) 2020 Richard Thompson, Qatar Biomedical Research Institute
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
 *
 */
 
 #include <seqan3/argument_parser/all.hpp> // for argument_parser
#include <seqan3/core/debug_stream.hpp>   // for debug_stream
#include <seqan3/std/filesystem>                          // for tmp_dir
#include <seqan3/io/sequence_file/input.hpp>    // for sequence_file_input
 
int main()
{
    seqan3::debug_stream << "Hello world\n";
    return 0;
}
