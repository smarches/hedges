#include <cstddef>
#include <string>
#include <memory>
#include <iostream>
#include <vector>
#include <array>

#include "schifra/schifra_galois_field.hpp"
#include "schifra/schifra_galois_field_polynomial.hpp"
#include "schifra/schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra/schifra_reed_solomon_encoder.hpp"
#include "schifra/schifra_reed_solomon_decoder.hpp"
#include "schifra/schifra_reed_solomon_block.hpp"
#include "schifra/schifra_error_processes.hpp"

using Uchar = std::uint8_t;
using VecUchar = std::vector<Uchar>;

union b8 {
	unsigned long long w;
	unsigned char b[8];
	b8() {}
	b8(long long ww) : w(ww) {}
	b8(const VecUchar &vec) {
		if (vec.size() != 8) throw("bad conversion of VecUchar to b8");
		for (int i = 0; i < 8; i++) b[i] = vec[i];
	}
	bool operator==(const b8 &b) { return w == b.w; }
	bool operator!=(const b8 &b) { return w != b.w; }
	bool operator>(const b8 &b)  { return w > b.w; }
	bool operator<(const b8 &b)  { return w < b.w; }
	bool operator>=(const b8 &b) { return w >= b.w; }
	bool operator<=(const b8 &b) { return w <= b.w; }
};

using namespace schifra;

template<size_t code_length, size_t fec_length>
class SchifraCode { // see Schifra's example08
    // see schifra_galois_field.hpp
    static_assert(code_length > fec_length);
	static constexpr size_t data_length = code_length - fec_length;
	static constexpr size_t generator_polynomial_root_count = fec_length;
	static constexpr int field_descriptor = 8;
	static constexpr size_t generator_polynomial_index = 120;
    static constexpr auto pad8 = std::min(8ul,data_length);
    static constexpr std::array<unsigned,9> pp6{ 1, 1, 1, 0, 0, 0, 0, 1, 1 };

	using encoder_t = reed_solomon::shortened_encoder<code_length, fec_length, data_length>;
	using decoder_t = reed_solomon::shortened_decoder<code_length, fec_length, data_length>;
	using field_t = galois::field;
	using field_polynomial_t = galois::field_polynomial;
	using block_t = reed_solomon::block<code_length, fec_length>;

    // for some reason the default constructors in Schifra are
    // private so we need to declare pointers to the things that will be
    // initialized, and all the subsequent code is poisoned
    field_t *field;
	encoder_t *encoder;
	decoder_t *decoder;
	field_polynomial_t *generator_polynomial;
	block_t *block; // size of block-> data is template parameter block_length
    // return value info of decode()
    template<typename T>
    struct rv {
		int err_number, errs_detected, errs_corrected;
        bool recoverable;
		T message;
        rv() = delete;
        rv(block_t* b) : err_number(b->error),
            errs_detected(b->errors_detected),
		    errs_corrected(b->errors_corrected),
		    recoverable(!b->unrecoverable) {}
    };
    
    void _encode(std::string& data) {
        auto res = encoder->encode(data, *block);
        if (!res) {
            std::cerr << "Error - Critical encoding failure!\nMsg: " << block->error_as_string() << '\n';
            exit(EXIT_FAILURE);
        }
    }

public:

	SchifraCode() {
        /* Instantiate Finite Field, Generator Polynomials, Encoder, Decoder */
        // field = new field_t(field_descriptor,pp6.size(), pp6.data());
        field = new field_t(field_descriptor,pp6.size(), pp6.data());
		generator_polynomial = new field_polynomial_t(*field);
        auto res = make_sequential_root_generator_polynomial(
            *field,
            generator_polynomial_index,
            generator_polynomial_root_count,
            *generator_polynomial
        );
		if ( !res ) { throw("Error - Failed to create sequential root generator!"); }
		encoder = new encoder_t(*field, *generator_polynomial);
		decoder = new decoder_t(*field, generator_polynomial_index);
	}

    ~SchifraCode() {
        delete field;
        delete generator_polynomial;
        delete encoder;
        delete decoder;
    }

	VecUchar encode(VecUchar& msg) {
		std::string message(code_length, 0);
        std::size_t L = std::min(msg.size(),data_length);
		for (size_t i = 0; i < L; i++) message[i] = msg[i];
        _encode(message);
		VecUchar codeword(code_length);
		for (int i = 0; i < code_length; i++) codeword[i] = static_cast<Uchar>(block->data[i]);
		return codeword;
	}
	b8 encode(b8 msg) {
		std::string message(data_length,0);
        // 'msg.b' is exactly 8 bytes
		for (unsigned i = 0; i < 8; i++) message.push_back(msg.b[i]);
        // the encode_t.encode() method blindly copies the input message up to its data_length
        // template parameter regardless of whether the message is actually that long
		b8 codeword(0);
		for (size_t i = 0; i < pad8; i++) codeword.b[i] = static_cast<Uchar>(block->data[i]);
		return codeword;
	}
    
	rv<VecUchar> decode(VecUchar& codeword) {
        if(codeword.size() != code_length) {
            throw("cannot decode this!");
        }
		block->reset(); // not needed?
		for (int i = 0; i < code_length; i++) {
            block->data[i] = static_cast<galois::field_symbol>(codeword[i]);
        }
		decoder.decode(block);
		VecUchar message(code_length);
		for (int i = 0; i < code_length; i++) message[i] = static_cast<Uchar>(block->data[i]);
        rv<VecUchar> ret(block);
        ret.message(std::move(message));
		return ret;
	}
	rv<VecUchar> decode(const VecUchar& codeword,const std::vector<size_t>& erasures) {
        if(codeword.size() != code_length) {
            throw("cannot decode this!");
        }
		block->reset(); // not needed?
		for (int i = 0; i < code_length; i++) {
            block->data[i] = static_cast<galois::field_symbol>(codeword[i]);
        }
		decoder->decode(block,erasures);
		rv<VecUchar> ret(block);
        VecUchar message(code_length);
		for (int i = 0; i < code_length; i++) message[i] = static_cast<Uchar>(block->data[i]);
        ret.message(std::move(message));
		return ret;
	}
	auto decode(b8 codeword) {
		block->reset(); // not needed?
		for (size_t i = 0; i < code_length; i++) {
            block->data[i] = static_cast<galois::field_symbol>(codeword.b[i]);
        }
		decoder->decode(*block);
        std::cout << "actual decode finished\n";
		b8 message(*(unsigned long long *)(block->data));
        std::cout << "b8 message is: " << message.w << '\n';
        rv<b8> ret(block);
        ret.message = message;
		return ret;
	}
	auto decode(const b8 codeword,const std::vector<size_t>& erasures) {
		block->reset(); // not needed?
		for (int i = 0; i < code_length; i++){
            block->data[i] = static_cast<galois::field_symbol>(codeword.b[i]);
        }
		decoder->decode(*block,erasures);
		b8 message(*(unsigned long long *)(block->data));
        std::cout << "b8 message is: " << message.w << '\n';
        rv<b8> ret(block);
        ret.message = message;
		return ret;
	}

	/*
	0 "No Error";
	1 "Invalid Encoder";
	2 "Incompatible Generator Polynomial";
	3 "Invalid Decoder";
	4 "Decoder Failure - Non-zero Syndrome";
	5 "Decoder Failure - Too Many Errors/Erasures";
	6 "Decoder Failure - Invalid Symbol Correction";
	7 "Decoder Failure - Invalid Codeword Correction";
	*/
	size_t code_len() const noexcept { return code_length; }
	size_t fec_len()  const noexcept { return fec_length; }
	size_t data_len() const noexcept { return data_length; }
};
