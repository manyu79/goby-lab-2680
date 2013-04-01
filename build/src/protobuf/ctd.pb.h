// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: ctd.proto

#ifndef PROTOBUF_ctd_2eproto__INCLUDED
#define PROTOBUF_ctd_2eproto__INCLUDED

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
void  protobuf_AddDesc_ctd_2eproto();
void protobuf_AssignDesc_ctd_2eproto();
void protobuf_ShutdownFile_ctd_2eproto();

class CTDMessage;

// ===================================================================

class CTDMessage : public ::google::protobuf::Message {
 public:
  CTDMessage();
  virtual ~CTDMessage();
  
  CTDMessage(const CTDMessage& from);
  
  inline CTDMessage& operator=(const CTDMessage& from) {
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
  static const CTDMessage& default_instance();
  
  void Swap(CTDMessage* other);
  
  // implements Message ----------------------------------------------
  
  CTDMessage* New() const;
  void CopyFrom(const ::google::protobuf::Message& from);
  void MergeFrom(const ::google::protobuf::Message& from);
  void CopyFrom(const CTDMessage& from);
  void MergeFrom(const CTDMessage& from);
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
  
  // accessors -------------------------------------------------------
  
  // repeated int32 depth = 1;
  inline int depth_size() const;
  inline void clear_depth();
  static const int kDepthFieldNumber = 1;
  inline ::google::protobuf::int32 depth(int index) const;
  inline void set_depth(int index, ::google::protobuf::int32 value);
  inline void add_depth(::google::protobuf::int32 value);
  inline const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
      depth() const;
  inline ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
      mutable_depth();
  
  // repeated int32 temperature = 2;
  inline int temperature_size() const;
  inline void clear_temperature();
  static const int kTemperatureFieldNumber = 2;
  inline ::google::protobuf::int32 temperature(int index) const;
  inline void set_temperature(int index, ::google::protobuf::int32 value);
  inline void add_temperature(::google::protobuf::int32 value);
  inline const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
      temperature() const;
  inline ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
      mutable_temperature();
  
  // repeated double salinity = 3;
  inline int salinity_size() const;
  inline void clear_salinity();
  static const int kSalinityFieldNumber = 3;
  inline double salinity(int index) const;
  inline void set_salinity(int index, double value);
  inline void add_salinity(double value);
  inline const ::google::protobuf::RepeatedField< double >&
      salinity() const;
  inline ::google::protobuf::RepeatedField< double >*
      mutable_salinity();
  
  // @@protoc_insertion_point(class_scope:CTDMessage)
 private:
  
  ::google::protobuf::UnknownFieldSet _unknown_fields_;
  
  ::google::protobuf::RepeatedField< ::google::protobuf::int32 > depth_;
  ::google::protobuf::RepeatedField< ::google::protobuf::int32 > temperature_;
  ::google::protobuf::RepeatedField< double > salinity_;
  
  mutable int _cached_size_;
  ::google::protobuf::uint32 _has_bits_[(3 + 31) / 32];
  
  friend void  protobuf_AddDesc_ctd_2eproto();
  friend void protobuf_AssignDesc_ctd_2eproto();
  friend void protobuf_ShutdownFile_ctd_2eproto();
  
  void InitAsDefaultInstance();
  static CTDMessage* default_instance_;
};
// ===================================================================


// ===================================================================

// CTDMessage

// repeated int32 depth = 1;
inline int CTDMessage::depth_size() const {
  return depth_.size();
}
inline void CTDMessage::clear_depth() {
  depth_.Clear();
}
inline ::google::protobuf::int32 CTDMessage::depth(int index) const {
  return depth_.Get(index);
}
inline void CTDMessage::set_depth(int index, ::google::protobuf::int32 value) {
  depth_.Set(index, value);
}
inline void CTDMessage::add_depth(::google::protobuf::int32 value) {
  depth_.Add(value);
}
inline const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
CTDMessage::depth() const {
  return depth_;
}
inline ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
CTDMessage::mutable_depth() {
  return &depth_;
}

// repeated int32 temperature = 2;
inline int CTDMessage::temperature_size() const {
  return temperature_.size();
}
inline void CTDMessage::clear_temperature() {
  temperature_.Clear();
}
inline ::google::protobuf::int32 CTDMessage::temperature(int index) const {
  return temperature_.Get(index);
}
inline void CTDMessage::set_temperature(int index, ::google::protobuf::int32 value) {
  temperature_.Set(index, value);
}
inline void CTDMessage::add_temperature(::google::protobuf::int32 value) {
  temperature_.Add(value);
}
inline const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
CTDMessage::temperature() const {
  return temperature_;
}
inline ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
CTDMessage::mutable_temperature() {
  return &temperature_;
}

// repeated double salinity = 3;
inline int CTDMessage::salinity_size() const {
  return salinity_.size();
}
inline void CTDMessage::clear_salinity() {
  salinity_.Clear();
}
inline double CTDMessage::salinity(int index) const {
  return salinity_.Get(index);
}
inline void CTDMessage::set_salinity(int index, double value) {
  salinity_.Set(index, value);
}
inline void CTDMessage::add_salinity(double value) {
  salinity_.Add(value);
}
inline const ::google::protobuf::RepeatedField< double >&
CTDMessage::salinity() const {
  return salinity_;
}
inline ::google::protobuf::RepeatedField< double >*
CTDMessage::mutable_salinity() {
  return &salinity_;
}


// @@protoc_insertion_point(namespace_scope)

#ifndef SWIG
namespace google {
namespace protobuf {


}  // namespace google
}  // namespace protobuf
#endif  // SWIG

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_ctd_2eproto__INCLUDED
