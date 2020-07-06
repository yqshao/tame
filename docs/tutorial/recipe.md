# Writing a recipe

## Recipes in TAME

A recipe is essentially a python program that runs specific tasks.
Below is how a typical recipe looks like:

```Python
   def actual_fn(arg1, arg2=1, ...):
       do some actual work

   def set_parser(parser):
       parser.add_argument('arg1', help='argument 1')
       parser.add_argument('--arg2', help='argument 2', default=1)
       parser.set_defaults(func=lambda args: actual_fn(
           args.arg1, args.arg2))

   def main():
       import argparse
       parser = argparse.ArgumentParser()
       set_parser(parser)
       args = parser.parse_args()
       args.func(args)

   if __name__ == "__main__":
       main()
```

The recipes will be gathered as subcommands of `tame`. To do so,
each tame recipes needs a `set_parser` function which sets the
command line arguments and documentations (see one of the recipes for
an idea). `set_parser` will be called both when the reciepe is used
as a standalone script or as a subcommnad.

## Checklist before contributing a recipe

Here is an additional checklist for implementing a new recipe.

- place the recipe in a sub-package of recipes.
- test the recipe with the latest tame code
- document the method used and cite appropriate reference in docstring
  of the recipe.