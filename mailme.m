function mailme(~,mailtitle,mailcontent)
% 账号设置
mail = '1900011012@pku.edu.cn';  % ①发送邮件的邮箱地址
password = 'zspzhousp123'; % ②发送邮件邮箱授权码
% 服务器设置
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.pku.edu.cn'); % ③SMTP服务器，这里我选择了QQ邮箱
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
% 发送邮件
receiver='1772592317@qq.com'; % ④我的收件邮箱。可以设为缺省值或不设
sendmail(receiver,mailtitle,mailcontent);
end