from django import template
from django.template.defaultfilters import stringfilter
from django.template.loader import get_template
import re

register = template.Library()

@register.filter
@stringfilter
def split_text(value,length=60):
	wrapped = ''
	length = int(length)
	for x in range(0,len(value),length):
		wrapped+=str(value[x:x+length])+'\n'

	wrapped = wrapped.rstrip()
	return wrapped


@register.simple_tag
def help_tip(id):
        
    return '<a tabindex="0" role="button" data-toggle="popover"'\
    		+' data-trigger="focus" id="'+id+'_trigger"'\
    		+'><i class="glyphicon glyphicon-question-sign"></i></a>';


