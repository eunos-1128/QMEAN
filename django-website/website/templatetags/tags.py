from django import template
from django.template.defaultfilters import stringfilter

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