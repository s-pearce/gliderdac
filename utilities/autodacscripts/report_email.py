import configparser as cp
import sys
import pdb
import smtplib
import datetime as dt
import email
from email.mime.text import MIMEText

class ReportEmail(object):
    """A class to report error messages via email to recipients from the 
    configurations file from the Glider DAC code.
    """
    def __init__(self, filename=None):
        self.recipients = []
        self.config = None
        self.hostname = ''
        self.sender = ''
        self.port = 25
        self.subject = 'DAC FTP upload error {:s}'
        #self.get_config_info(filename)

    def get_config_info(self, filename):
        self.config = cp.RawConfigParser()
        with open(filename, 'r') as configfile:
            self.config.readfp(configfile)
        recipients = self.config.get(section, field)
        if recipients:
            self.recipients = recipients.split(', ')
        sender = self.config.get(section, field)
        if sender:
            self.sender = sender
        hostname = self.config.get(section, field)
        if hostname:
            self.hostname = hostname
        port = self.config.get(section, field)
        if port:
            self.port = port

    def sendmail(self, msg):
        success = None
        now = dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        subject = self.subject.format(now)
        msg_prefix = "From:{:s}\n".format(self.sender)
        msg_prefix += "To:{:s}\n".format(','.join(self.recipients))
        msg_prefix += "Subject:{:s}\n\n".format(subject)
        msg = msg_prefix + msg

        try:
            smtp_server = smtplib.SMTP(host=self.hostname, port=self.port)
            s = smtp_server.sendmail(self.sender, self.recipients, msg)
        except smtplib.SMTPException:
            success = False
        else:
            success = True
        finally:
            smtp_server.quit()
        return (success, s)