// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: acomms_example.proto

#ifndef PROTOBUF_acomms_5fexample_2eproto__INCLUDED
#define PROTOBUF_acomms_5fexample_2eproto__INCLUDED

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 2004000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 2004001 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/repeated_field.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/generated_message_reflection.h>
#include "goby/acomms/protobuf/dccl_option_extensions.pb.h"
// @@protoc_insertion_point(includes)

// Internal implementation detail -- do not call these.
void  protobuf_AddDesc_acomms_5fexample_2eproto();
void protobuf_AssignDesc_acomms_5fexample_2eproto();
void protobuf_ShutdownFile_acomms_5fexample_2eproto();

class AcommsExampleMessage;

enum AcommsExampleMessage_Sender {
  AcommsExampleMessage_Sender_HENRIK = 1,
  AcommsExampleMessage_Sender_MIKERB = 2,
  AcommsExampleMessage_Sender_TES = 3,
  AcommsExampleMessage_Sender_AYAARI = 4,
  AcommsExampleMessage_Sender_ARTHURA = 5,
  AcommsExampleMessage_Sender_TBARTON = 6,
  AcommsExampleMessage_Sender_CLOITRE = 7,
  AcommsExampleMessage_Sender_SDANESH = 8,
  AcommsExampleMessage_Sender_INE = 9,
  AcommsExampleMessage_Sender_EMF43 = 10,
  AcommsExampleMessage_Sender_ARG = 11,
  AcommsExampleMessage_Sender_MATTGMIT = 12,
  AcommsExampleMessage_Sender_JLEIGHT = 13,
  AcommsExampleMessage_Sender_MMCGRAW = 14,
  AcommsExampleMessage_Sender_SPETILLO = 15,
  AcommsExampleMessage_Sender_MANGPO = 16,
  AcommsExampleMessage_Sender_ICRUST = 17,
  AcommsExampleMessage_Sender_ABS198 = 18,
  AcommsExampleMessage_Sender_JSCHUL = 19,
  AcommsExampleMessage_Sender_ROBTRUAX = 20,
  AcommsExampleMessage_Sender_YCK = 21
};
bool AcommsExampleMessage_Sender_IsValid(int value);
const AcommsExampleMessage_Sender AcommsExampleMessage_Sender_Sender_MIN = AcommsExampleMessage_Sender_HENRIK;
const AcommsExampleMessage_Sender AcommsExampleMessage_Sender_Sender_MAX = AcommsExampleMessage_Sender_YCK;
const int AcommsExampleMessage_Sender_Sender_ARRAYSIZE = AcommsExampleMessage_Sender_Sender_MAX + 1;

const ::google::protobuf::EnumDescriptor* AcommsExampleMessage_Sender_descriptor();
inline const ::std::string& AcommsExampleMessage_Sender_Name(AcommsExampleMessage_Sender value) {
  return ::google::protobuf::internal::NameOfEnum(
    AcommsExampleMessage_Sender_descriptor(), value);
}
inline bool AcommsExampleMessage_Sender_Parse(
    const ::std::string& name, AcommsExampleMessage_Sender* value) {
  return ::google::protobuf::internal::ParseNamedEnum<AcommsExampleMessage_Sender>(
    AcommsExampleMessage_Sender_descriptor(), name, value);
}
// ===================================================================

class AcommsExampleMessage : public ::google::protobuf::Message {
 public:
  AcommsExampleMessage();
  virtual ~AcommsExampleMessage();
  
  AcommsExampleMessage(const AcommsExampleMessage& from);
  
  inline AcommsExampleMessage& operator=(const AcommsExampleMessage& from) {
    CopyFrom(from);
    return *this;
  }
  
  inline const ::google::protobuf::UnknownFieldSet& unknown_fields() const {
    return _unknown_fields_;
  }
  
  inline ::google::protobuf::UnknownFieldSet* mutable_unknown_fields() {
    return &_unknown_fields_;
  }
  
  static const ::google::protobuf::Descriptor* descriptor();
  static const AcommsExampleMessage& default_instance();
  
  void Swap(AcommsExampleMessage* other);
  
  // implements Message ----------------------------------------------
  
  AcommsExampleMessage* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const AcommsExampleMessage& from);
  void MergeFrom(const AcommsExampleMessage& from);
  void Clear();
  bool IsInitialized() const;
  
  int ByteSize() const;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input);
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const;
  ::google::protobuf::uint8* SerializeWithCachedSizesToArray(::google::protobuf::uint8* output) const;
  int GetCachedSize() const { return _cached_size_; }
  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const;
  public:
  
  ::google::protobuf::Metadata GetMetadata() const;
  
  // nested types ----------------------------------------------------
  
  typedef AcommsExampleMessage_Sender Sender;
  static const Sender HENRIK = AcommsExampleMessage_Sender_HENRIK;
  static const Sender MIKERB = AcommsExampleMessage_Sender_MIKERB;
  static const Sender TES = AcommsExampleMessage_Sender_TES;
  static const Sender AYAARI = AcommsExampleMessage_Sender_AYAARI;
  static const Sender ARTHURA = AcommsExampleMessage_Sender_ARTHURA;
  static const Sender TBARTON = AcommsExampleMessage_Sender_TBARTON;
  static const Sender CLOITRE = AcommsExampleMessage_Sender_CLOITRE;
  static const Sender SDANESH = AcommsExampleMessage_Sender_SDANESH;
  static const Sender INE = AcommsExampleMessage_Sender_INE;
  static const Sender EMF43 = AcommsExampleMessage_Sender_EMF43;
  static const Sender ARG = AcommsExampleMessage_Sender_ARG;
  static const Sender MATTGMIT = AcommsExampleMessage_Sender_MATTGMIT;
  static const Sender JLEIGHT = AcommsExampleMessage_Sender_JLEIGHT;
  static const Sender MMCGRAW = AcommsExampleMessage_Sender_MMCGRAW;
  static const Sender SPETILLO = AcommsExampleMessage_Sender_SPETILLO;
  static const Sender MANGPO = AcommsExampleMessage_Sender_MANGPO;
  static const Sender ICRUST = AcommsExampleMessage_Sender_ICRUST;
  static const Sender ABS198 = AcommsExampleMessage_Sender_ABS198;
  static const Sender JSCHUL = AcommsExampleMessage_Sender_JSCHUL;
  static const Sender ROBTRUAX = AcommsExampleMessage_Sender_ROBTRUAX;
  static const Sender YCK = AcommsExampleMessage_Sender_YCK;
  static inline bool Sender_IsValid(int value) {
    return AcommsExampleMessage_Sender_IsValid(value);
  }
  static const Sender Sender_MIN =
    AcommsExampleMessage_Sender_Sender_MIN;
  static const Sender Sender_MAX =
    AcommsExampleMessage_Sender_Sender_MAX;
  static const int Sender_ARRAYSIZE =
    AcommsExampleMessage_Sender_Sender_ARRAYSIZE;
  static inline const ::google::protobuf::EnumDescriptor*
  Sender_descriptor() {
    return AcommsExampleMessage_Sender_descriptor();
  }
  static inline const ::std::string& Sender_Name(Sender value) {
    return AcommsExampleMessage_Sender_Name(value);
  }
  static inline bool Sender_Parse(const ::std::string& name,
      Sender* value) {
    return AcommsExampleMessage_Sender_Parse(name, value);
  }
  
  // accessors -------------------------------------------------------
  
  // required uint32 index = 1;
  inline bool has_index() const;
  inline void clear_index();
  static const int kIndexFieldNumber = 1;
  inline ::google::protobuf::uint32 index() const;
  inline void set_index(::google::protobuf::uint32 value);
  
  // required int32 random = 2;
  inline bool has_random() const;
  inline void clear_random();
  static const int kRandomFieldNumber = 2;
  inline ::google::protobuf::int32 random() const;
  inline void set_random(::google::protobuf::int32 value);
  
  // required .AcommsExampleMessage.Sender sender = 3;
  inline bool has_sender() const;
  inline void clear_sender();
  static const int kSenderFieldNumber = 3;
  inline ::AcommsExampleMessage_Sender sender() const;
  inline void set_sender(::AcommsExampleMessage_Sender value);
  
  // @@protoc_insertion_point(class_scope:AcommsExampleMessage)
 private:
  inline void set_has_index();
  inline void clear_has_index();
  inline void set_has_random();
  inline void clear_has_random();
  inline void set_has_sender();
  inline void clear_has_sender();
  
  ::google::protobuf::UnknownFieldSet _unknown_fields_;
  
  ::google::protobuf::uint32 index_;
  ::google::protobuf::int32 random_;
  int sender_;
  
  mutable int _cached_size_;
  ::google::protobuf::uint32 _has_bits_[(3 + 31) / 32];
  
  friend void  protobuf_AddDesc_acomms_5fexample_2eproto();
  friend void protobuf_AssignDesc_acomms_5fexample_2eproto();
  friend void protobuf_ShutdownFile_acomms_5fexample_2eproto();
  
  void InitAsDefaultInstance();
  static AcommsExampleMessage* default_instance_;
};
// ===================================================================


// ===================================================================

// AcommsExampleMessage

// required uint32 index = 1;
inline bool AcommsExampleMessage::has_index() const {
  return (_has_bits_[0] & 0x00000001u) != 0;
}
inline void AcommsExampleMessage::set_has_index() {
  _has_bits_[0] |= 0x00000001u;
}
inline void AcommsExampleMessage::clear_has_index() {
  _has_bits_[0] &= ~0x00000001u;
}
inline void AcommsExampleMessage::clear_index() {
  index_ = 0u;
  clear_has_index();
}
inline ::google::protobuf::uint32 AcommsExampleMessage::index() const {
  return index_;
}
inline void AcommsExampleMessage::set_index(::google::protobuf::uint32 value) {
  set_has_index();
  index_ = value;
}

// required int32 random = 2;
inline bool AcommsExampleMessage::has_random() const {
  return (_has_bits_[0] & 0x00000002u) != 0;
}
inline void AcommsExampleMessage::set_has_random() {
  _has_bits_[0] |= 0x00000002u;
}
inline void AcommsExampleMessage::clear_has_random() {
  _has_bits_[0] &= ~0x00000002u;
}
inline void AcommsExampleMessage::clear_random() {
  random_ = 0;
  clear_has_random();
}
inline ::google::protobuf::int32 AcommsExampleMessage::random() const {
  return random_;
}
inline void AcommsExampleMessage::set_random(::google::protobuf::int32 value) {
  set_has_random();
  random_ = value;
}

// required .AcommsExampleMessage.Sender sender = 3;
inline bool AcommsExampleMessage::has_sender() const {
  return (_has_bits_[0] & 0x00000004u) != 0;
}
inline void AcommsExampleMessage::set_has_sender() {
  _has_bits_[0] |= 0x00000004u;
}
inline void AcommsExampleMessage::clear_has_sender() {
  _has_bits_[0] &= ~0x00000004u;
}
inline void AcommsExampleMessage::clear_sender() {
  sender_ = 1;
  clear_has_sender();
}
inline ::AcommsExampleMessage_Sender AcommsExampleMessage::sender() const {
  return static_cast< ::AcommsExampleMessage_Sender >(sender_);
}
inline void AcommsExampleMessage::set_sender(::AcommsExampleMessage_Sender value) {
  GOOGLE_DCHECK(::AcommsExampleMessage_Sender_IsValid(value));
  set_has_sender();
  sender_ = value;
}


// @@protoc_insertion_point(namespace_scope)

#ifndef SWIG
namespace google {
namespace protobuf {

template <>
inline const EnumDescriptor* GetEnumDescriptor< ::AcommsExampleMessage_Sender>() {
  return ::AcommsExampleMessage_Sender_descriptor();
}

}  // namespace google
}  // namespace protobuf
#endif  // SWIG

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_acomms_5fexample_2eproto__INCLUDED