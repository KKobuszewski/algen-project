#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// xml parsing
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

// if libxml2 not usable will throw error during the compilation due to undefined macro ENCODING
#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)
#define ENCODING "UTF-8"

#define XML_BUFFER_SIZE 256
#define MAX_STRING_SIZE 1024

namespace xml_raport
{

int rc;
xmlTextWriterPtr writer;
xmlChar tmp_names_in[XML_BUFFER_SIZE][MAX_STRING_SIZE];
xmlChar tmp_data_in[XML_BUFFER_SIZE][MAX_STRING_SIZE];
xmlChar tmp_names_out[XML_BUFFER_SIZE][MAX_STRING_SIZE];
xmlChar tmp_data_out[XML_BUFFER_SIZE][MAX_STRING_SIZE];
unsigned it_in = 0;
unsigned it_out = 0;


inline void write_input(const char* name, const double datum)
{
    if (it_in >= XML_BUFFER_SIZE)  { fprintf(stderr,"Error! XML input buffer overloaded!\n"); exit(EXIT_FAILURE); }
    
    strcpy( reinterpret_cast< char *>(tmp_names_in[it_in]),name );
    if ( snprintf( reinterpret_cast< char *>(tmp_data_in[it_in]),MAX_STRING_SIZE,"%.15e",datum ) > MAX_STRING_SIZE )
    { fprintf(stderr,"Error! Char buffer overloaded!"); exit(EXIT_FAILURE); }
    
    it_in++; //
}

inline void write_output(const char* name, const double datum)
{
    if (it_out >= XML_BUFFER_SIZE)  { fprintf(stderr,"Error! XML output buffer overloaded!\n"); exit(EXIT_FAILURE); }
    
    strcpy( reinterpret_cast< char *>(tmp_names_out[it_out]),name );
    if ( snprintf( reinterpret_cast< char *>(tmp_data_out[it_out]),MAX_STRING_SIZE,"%.15e",datum ) > MAX_STRING_SIZE )
    { fprintf(stderr,"Error! Char buffer overloaded!"); exit(EXIT_FAILURE); }
    
    it_out++; //
}

inline void save_raport(const char* uri)
{
    // ======================================= Create xml file ===========================================
    /* Create a new XmlWriter for uri, with no compression. */
    writer = xmlNewTextWriterFilename(uri, 0);
    if (writer == NULL) {
        printf("testXmlwriterFilename: Error creating the xml writer\n");
        exit(EXIT_FAILURE);
    }
    
    /* Start the document with the xml default for the version,
     * encoding UTF-8 and the default for the standalone
     * declaration. */
    rc = xmlTextWriterStartDocument(writer, NULL, ENCODING, NULL);
    if (rc < 0) {
        printf
            ("testXmlwriterFilename: Error at xmlTextWriterStartDocument\n");
        exit(EXIT_FAILURE);
    }
    
    // ======================================= Create root ===============================================
    /* Start an element named "EXAMPLE". Since thist is the first
     * element, this will be the root element of the document. */
    rc = xmlTextWriterStartElement(writer, BAD_CAST "DATASHEET");
    if (rc < 0) {
        printf
            ("testXmlwriterFilename: Error at xmlTextWriterStartElement\n");
        exit(EXIT_FAILURE);
    }
    
    
    // ==================================== Write input to xml ===============================================
    /* Start an element named "INPUT". */
    rc = xmlTextWriterStartElement(writer, BAD_CAST "INPUT");
    if (rc < 0) {
        printf
            ("testXmlwriterFilename: Error at xmlTextWriterStartElement\n");
        exit(EXIT_FAILURE);
    }
    
    
    for (unsigned ii=0; ii < it_in; ii++)
    {
        rc = xmlTextWriterStartElement(writer, BAD_CAST "ENTRY");
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterStartElement\n");
            exit(EXIT_FAILURE);
        }
        
        /* Write an element named "NAME" as child of ENTRY. */
        rc = xmlTextWriterWriteElement(writer, BAD_CAST "NAME",
                                    BAD_CAST tmp_names_in[ii]);
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterWriteElement\n");
            return;
        }
        /* Write an element named "NAME" as child of ENTRY. */
        rc = xmlTextWriterWriteElement(writer, BAD_CAST "VALUE",
                                    BAD_CAST tmp_data_in[ii]);
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterWriteElement\n");
            return;
        }
        
        
        /* Close the element named ENTRY. */
        rc = xmlTextWriterEndElement(writer);
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterEndElement\n");
            exit(EXIT_FAILURE);
        }
    }
    /* Close the element named INPUT. */
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf
            ("testXmlwriterFilename: Error at xmlTextWriterEndElement\n");
        exit(EXIT_FAILURE);
    }
    
    
    
    // ======================================= Write output to xml ============================================
    /* Start an element named "OUTPUT". */
    rc = xmlTextWriterStartElement(writer, BAD_CAST "OUTPUT");
    if (rc < 0) {
        printf("testXmlwriterFilename: Error at xmlTextWriterStartElement\n");
        exit(EXIT_FAILURE);
    }
    
    
    for (unsigned ii=0; ii < it_out; ii++)
    {
        rc = xmlTextWriterStartElement(writer, BAD_CAST "ENTRY");
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterStartElement\n");
            exit(EXIT_FAILURE);
        }
        
        /* Write an element named "NAME" as child of ENTRY. */
        rc = xmlTextWriterWriteElement(writer, BAD_CAST "NAME",
                                    BAD_CAST tmp_names_out[ii]);
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterWriteElement\n");
            return;
        }
        /* Write an element named "NAME" as child of ENTRY. */
        rc = xmlTextWriterWriteElement(writer, BAD_CAST "VALUE",
                                    BAD_CAST tmp_data_out[ii]);
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterWriteElement\n");
            return;
        }
        
        
        /* Close the element named ENTRY. */
        rc = xmlTextWriterEndElement(writer);
        if (rc < 0) {
            printf
                ("testXmlwriterFilename: Error at xmlTextWriterEndElement\n");
            exit(EXIT_FAILURE);
        }
    }
    /* Close the element named OUTPUT. */
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf
            ("testXmlwriterFilename: Error at xmlTextWriterEndElement\n");
        exit(EXIT_FAILURE);
    }
    
    
    
    
    // ================================== End document ================================================
     
     
    /* Close the element named DATASHEET. */
    rc = xmlTextWriterEndElement(writer);
    if (rc < 0) {
        printf
            ("testXmlwriterFilename: Error at xmlTextWriterEndElement\n");
        exit(EXIT_FAILURE);
    }
     
     
    /* Here we could close the elements using the
     * function xmlTextWriterEndElement, but since we do not want to
     * write any other elements, we simply call xmlTextWriterEndDocument,
     * which will do all the work. */
    rc = xmlTextWriterEndDocument(writer);
    if (rc < 0) {
        printf
            ("testXmlwriterFilename: Error at xmlTextWriterEndDocument\n");
        exit(EXIT_FAILURE);
    }

    xmlFreeTextWriter(writer);
}


}
#endif