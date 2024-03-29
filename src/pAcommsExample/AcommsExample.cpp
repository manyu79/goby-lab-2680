/************************************************************/
/*    NAME: Toby Schneider                                  */
/*    ORGN: MIT                                             */
/*    FILE: AcommsExample.cpp                               */
/*    DATE: 2012-03-12                                      */
/************************************************************/

#include <cstdlib> // for rand()
#include "goby/moos/moos_protobuf_helpers.h" // for goby::moos::MOOSTranslation
#include "goby/acomms/dccl.h" // for goby::acomms::DCCLCodec
#include "codecs/codec_load.h"

#include "AcommsExample.h"

//---------------------------------------------------------
// Constructor

AcommsExample::AcommsExample() :
    start_time_(MOOSTime()),
    post_interval_(2)
{
    goby::acomms::DCCLCodec* dccl = goby::acomms::DCCLCodec::get();
    goby_dccl_load(dccl);
    dccl->validate<AcommsExampleMessage>();
}

//---------------------------------------------------------
// Procedure: OnNewMail

bool AcommsExample::OnNewMail(MOOSMSG_LIST &NewMail)
{
   
    for(MOOSMSG_LIST::iterator it=NewMail.begin(), end=NewMail.end();
        it!=end;
        it++)
    {
        CMOOSMsg& msg = *it;
    }
	
    return(true);
}

//---------------------------------------------------------
// Procedure: OnConnectToServer

bool AcommsExample::OnConnectToServer()
{
    m_MissionReader.GetConfigurationParam("out_moos_var", out_moos_var_);	
    m_MissionReader.GetConfigurationParam("post_interval", post_interval_);	
    return(true);
}

//---------------------------------------------------------
// Procedure: Iterate()

bool AcommsExample::Iterate()
{
    // keep track of the message count
    static int i = 0;
    
    // is it time to send a message?
    if(MOOSTime() > (i*post_interval_ + start_time_))
    {
        // create a Protobuf message (AcommsExampleMessage autogenerated from acomms_message.proto)
        AcommsExampleMessage msg;
        msg.set_index(i);
        msg.set_random(rand() % 100);
        msg.set_sender(AcommsExampleMessage::TES);


        // serialize the protobuf message into a string suitable for MOOS
        std::string text_format_serialized;
        using goby::moos::protobuf::TranslatorEntry;
        goby::moos::MOOSTranslation<TranslatorEntry::TECHNIQUE_PROTOBUF_TEXT_FORMAT>::serialize(
            &text_format_serialized, msg);


        
        // post the serialized string
        m_Comms.Notify(out_moos_var_, text_format_serialized);
        
        std::cout << "Wrote to: [" << out_moos_var_ << "]: ["
                  << text_format_serialized << "]" << std::endl; 
        std::cout << "Sizes (bytes): TextFormat: " << text_format_serialized.size()
                  << ", Native Protobuf: " << msg.ByteSize()
                  << ", DCCL: " <<  goby::acomms::DCCLCodec::get()->size(msg) << std::endl;
        ++i;
    }
    

    return(true);
}

